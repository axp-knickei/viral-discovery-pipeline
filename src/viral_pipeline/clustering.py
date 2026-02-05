from __future__ import annotations

import subprocess
from pathlib import Path
from typing import List

import pandas as pd
import networkx as nx
from loguru import logger

from viral_pipeline.config import AppConfig
from viral_pipeline.utils import ClusteringResult, ensure_dir, check_dependency, group_paired_reads


class Clusterer:
    def __init__(self, config: AppConfig) -> None:
        self.config = config
        self.threads = str(config.tools.threads)

    def _run_command(self, cmd: List[str], description: str) -> None:
        logger.info(f"Running: {description}")
        logger.debug(f"Command: {' '.join(cmd)}")
        try:
            subprocess.run(
                cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
        except subprocess.CalledProcessError as e:
            logger.error(f"Error in {description}: {e.stderr}")
            raise RuntimeError(f"{description} failed") from e

    def cluster(self, contigs: Path, reads: List[Path], output_dir: Path) -> ClusteringResult:
        ensure_dir(output_dir)

        # Verify dependencies
        # tools: cd-hit-est, fastANI, bowtie2-build, bowtie2, samtools, coverm
        cd_hit = check_dependency("cd-hit-est")
        fast_ani = check_dependency("fastANI")
        bowtie2_build = check_dependency("bowtie2-build")
        bowtie2 = check_dependency("bowtie2")
        samtools = check_dependency("samtools")
        coverm = check_dependency("coverm")

        if not contigs.exists() or contigs.stat().st_size == 0:
            logger.warning("No viral contigs found. Skipping clustering and abundance.")
            # Return empty result
            clusters_file = output_dir / "clusters.tsv"
            abundance_file = output_dir / "vOTU_abundance_matrix_tpm.tsv"
            clusters_file.touch()
            abundance_file.touch()
            return ClusteringResult(clusters=clusters_file, abundance_table=abundance_file)

        # 1. Dereplication
        dereplicated_fa = output_dir / "all_viruses_dereplicated.fna"
        self._run_command(
            [
                str(cd_hit),
                "-i", str(contigs),
                "-o", str(dereplicated_fa),
                "-c", "0.99",
                "-n", "10",
                "-M", "0",
                "-T", self.threads
            ],
            "Dereplication with CD-HIT-EST"
        )

        # 2. Clustering (vOTUs)
        fastani_out = output_dir / "fastani_results.txt"
        self._run_command(
            [
                str(fast_ani),
                "--ql", str(dereplicated_fa),
                "--rl", str(dereplicated_fa),
                "-o", str(fastani_out),
                "-t", self.threads
            ],
            "Clustering calculation with FastANI"
        )

        # Parse FastANI and build graph
        votu_reps_fa = output_dir / "vOTUs_representatives.fa"
        clusters_tsv = output_dir / "clusters.tsv"
        self._process_clusters(
            dereplicated_fa, fastani_out, votu_reps_fa, clusters_tsv, samtools
        )

        # 3. Abundance Calculation
        abundance_matrix = self._calculate_abundance(
            votu_reps_fa, reads, output_dir, bowtie2_build, bowtie2, samtools, coverm
        )

        return ClusteringResult(clusters=clusters_tsv, abundance_table=abundance_matrix)

    def _process_clusters(
        self,
        dereplicated_fa: Path,
        fastani_out: Path,
        votu_reps_fa: Path,
        clusters_tsv: Path,
        samtools: Path
    ) -> None:
        logger.info("Processing clusters...")
        ani_threshold = self.config.clustering.ani_threshold
        min_coverage = self.config.clustering.min_coverage

        # Load FastANI results
        # Format: query, ref, ani, mapped_frags, total_frags
        try:
            df = pd.read_csv(
                fastani_out, sep="\t", header=None,
                names=["query", "ref", "ani", "mapped", "total"]
            )
        except pd.errors.EmptyDataError:
            logger.warning("FastANI produced no output. Treating all sequences as singletons.")
            df = pd.DataFrame(columns=["query", "ref", "ani", "mapped", "total"])

        # Filter edges
        # ANI is 0-100 in output usually? Bash script says $3/100 >= ani.
        # fastANI output is percentage (e.g., 99.5).
        edges = df[
            (df["ani"] / 100.0 >= ani_threshold) &
            (df["mapped"] / df["total"] >= min_coverage)
        ]

        # Build Graph
        G = nx.from_pandas_edgelist(edges, "query", "ref")

        # Add all nodes from the file (including singletons not in fastani output if any?)
        # FastANI usually outputs self-matches, so singletons should be there as query=ref.
        # But to be safe, we should ensure all sequences are in the graph.
        # We need to read the dereplicated FASTA headers.

        # Use samtools faidx to get sequence lengths
        self._run_command(
            [str(samtools), "faidx", str(dereplicated_fa)],
            "Indexing dereplicated FASTA"
        )

        seq_lengths = {}
        fai_file = Path(str(dereplicated_fa) + ".fai")
        with open(fai_file, "r") as f:
            for line in f:
                parts = line.split("\t")
                if len(parts) >= 2:
                    seq_lengths[parts[0]] = int(parts[1])

        # Efficiently add missing nodes (singletons)
        # add_nodes_from will ignore existing nodes if no attributes are provided,
        # preserving the graph structure while being faster than a Python loop.
        G.add_nodes_from(seq_lengths.keys())

        # Find components and select representatives (longest)
        reps = set()

        with open(clusters_tsv, "w") as f:
            f.write("cluster\tcontig\n")
            for i, component in enumerate(nx.connected_components(G)):
                # Sort by length descending
                sorted_nodes = sorted(component, key=lambda x: seq_lengths.get(x, 0), reverse=True)
                rep = sorted_nodes[0]
                reps.add(rep)
                votu_id = f"vOTU_{i+1:04d}" # Zero padded

                for node in sorted_nodes:
                    f.write(f"{votu_id}\t{node}\n")

        # Write representatives FASTA
        # Re-read file and write only reps
        with open(dereplicated_fa, "r") as fin, open(votu_reps_fa, "w") as fout:
            write_mode = False
            for line in fin:
                if line.startswith(">"):
                    seq_id = line[1:].split()[0]
                    if seq_id in reps:
                        write_mode = True
                        fout.write(line)
                    else:
                        write_mode = False
                elif write_mode:
                    fout.write(line)

        logger.info(f"Identified {len(reps)} vOTUs.")

    def _calculate_abundance(
        self,
        votu_db: Path,
        reads: List[Path],
        output_dir: Path,
        bowtie2_build: Path,
        bowtie2: Path,
        samtools: Path,
        coverm: Path
    ) -> Path:
        logger.info("Calculating abundance...")
        index_dir = output_dir / "index"
        ensure_dir(index_dir)
        index_base = index_dir / "vOTUs"

        # Build Index
        self._run_command(
            [str(bowtie2_build), str(votu_db), str(index_base)],
            "Building Bowtie2 index"
        )

        grouped_reads = group_paired_reads(reads)
        tmp_dir = output_dir / "tmp_abundance"
        ensure_dir(tmp_dir)

        tpm_files = []

        for sample_id, (r1, r2) in grouped_reads.items():
            bam_file = tmp_dir / f"{sample_id}.bam"
            tpm_file = tmp_dir / f"{sample_id}.tpm.txt"

            # Map
            cmd = [
                str(bowtie2),
                "-p", self.threads,
                "-x", str(index_base),
                "--very-sensitive"
            ]
            if r2:
                cmd.extend(["-1", str(r1), "-2", str(r2)])
            else:
                cmd.extend(["-U", str(r1)])

            # Pipe to samtools
            # Python subprocess piping: p1 | p2
            logger.info(f"Mapping reads for sample {sample_id}...")

            with open(bam_file, "wb") as bam_out:
                 # We pipe bowtie2 -> samtools view -> samtools sort -> file
                 # Or just bowtie2 | samtools sort -o file (samtools sort can read BAM from stdin if we give it -)
                 # Note: bowtie2 outputs SAM.

                 p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                 p2 = subprocess.Popen(
                     [str(samtools), "view", "-bS", "-"],
                     stdin=p1.stdout,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE
                 )
                 p3 = subprocess.Popen(
                     [str(samtools), "sort", "-@", self.threads, "-o", str(bam_file), "-"],
                     stdin=p2.stdout,
                     stdout=subprocess.PIPE,
                     stderr=subprocess.PIPE
                 )

                 out, err = p3.communicate()
                 p1.wait()
                 p2.wait()

                 if p3.returncode != 0:
                     logger.error(f"Mapping pipeline failed for {sample_id}: {err.decode()}")
                     raise RuntimeError(f"Mapping failed for {sample_id}")

            # CoverM
            logger.info(f"Running CoverM for {sample_id}...")
            # coverm contig -b input.bam -m tpm --min-covered-fraction 0.6
            self._run_command(
                [
                    str(coverm), "contig",
                    "-b", str(bam_file),
                    "-m", "tpm",
                    "--min-covered-fraction", "0.60",
                    "-o", str(tpm_file)
                ],
                f"CoverM for {sample_id}"
            )
            tpm_files.append(tpm_file)

        # Aggregate
        logger.info("Aggregating abundance matrix...")
        abundance_matrix_path = output_dir / "vOTU_abundance_matrix_tpm.tsv"

        dfs = []
        for f in tpm_files:
            # CoverM output: Genome\tTPM
            # We want to rename TPM column to sample ID
            sample_name = f.stem.replace(".tpm", "")
            d = pd.read_csv(f, sep="\t", header=0, index_col=0)
            d.rename(columns={d.columns[0]: sample_name}, inplace=True)
            dfs.append(d)

        if dfs:
            result_df = pd.concat(dfs, axis=1).fillna(0)
            result_df.index.name = "vOTU"
            result_df.to_csv(abundance_matrix_path, sep="\t")
        else:
            abundance_matrix_path.touch()

        return abundance_matrix_path
