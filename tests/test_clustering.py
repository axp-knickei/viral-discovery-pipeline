from unittest.mock import MagicMock
import pytest
from pathlib import Path
import pandas as pd
from viral_pipeline.clustering import Clusterer
from viral_pipeline.config import AppConfig

@pytest.fixture
def mock_config():
    # Helper to create config
    # We can rely on default values or mock objects
    config_dict = {
        "database": {
            "host_genome": "host.fa",
            "dvf_models": "dvf",
            "genomad_db": "genomad",
            "checkv_db": "checkv"
        },
        "tools": {
            "threads": 4,
            "memory_gb": 16,
            "min_contig_length": 1000
        },
        "quality": {
            "quality_threshold": 20,
            "min_read_length": 50,
            "adapter_removal": True
        },
        "clustering": {
            "ani_threshold": 0.95,
            "min_coverage": 0.85
        },
        "logging": {
            "level": "INFO",
            "log_dir": "logs"
        }
    }
    return AppConfig.model_validate(config_dict)

@pytest.fixture
def clusterer(mock_config):
    return Clusterer(mock_config)

def test_cluster_flow(clusterer, tmp_path, mocker):
    # Mocks
    mock_check_dependency = mocker.patch("viral_pipeline.clustering.check_dependency")
    mock_check_dependency.return_value = Path("/usr/bin/tool")

    mock_run = mocker.patch("viral_pipeline.clustering.subprocess.run")

    # Mock group_paired_reads to return a deterministic dict
    mocker.patch("viral_pipeline.clustering.group_paired_reads", return_value={
        "sample1": (Path("sample1_R1.fq"), Path("sample1_R2.fq"))
    })

    # Mock pandas read_csv for FastANI results
    # Query, Ref, ANI, Mapped, Total
    fastani_df = pd.DataFrame({
        "query": ["contig1", "contig1"],
        "ref": ["contig1", "contig2"],
        "ani": [100.0, 99.0],
        "mapped": [1000, 900],
        "total": [1000, 1000]
    })

    # We need to mock pd.read_csv to return fastani_df when called for fastani_results.txt
    # and to return coverm df when called for tpm.txt
    def mock_read_csv(filepath, *args, **kwargs):
        filepath = str(filepath)
        if "fastani_results.txt" in filepath:
            return fastani_df
        elif ".tpm.txt" in filepath:
            return pd.DataFrame({"Genome": ["vOTU_0001"], "TPM": [10.0]})
        return pd.DataFrame()

    mocker.patch("pandas.read_csv", side_effect=mock_read_csv)

    # Mock file inputs
    input_contigs = tmp_path / "contigs.fa"
    input_contigs.write_text(">contig1\nATCG\n>contig2\nATCG\n")

    input_reads = [Path("sample1_R1.fq"), Path("sample1_R2.fq")]
    output_dir = tmp_path / "output"

    # Mock subprocess Popen for the pipe chain (bowtie2 | samtools | samtools)
    mock_popen = mocker.patch("subprocess.Popen")
    process_mock = MagicMock()
    process_mock.communicate.return_value = (b"", b"")
    process_mock.returncode = 0
    mock_popen.return_value = process_mock

    # Simulate CD-HIT output (dereplicated file)
    ensure_dir = Path(output_dir)
    ensure_dir.mkdir(parents=True, exist_ok=True)
    derep_file = output_dir / "all_viruses_dereplicated.fna"
    derep_file.write_text(">contig1\nATCG\n>contig2\nATCG\n")

    # Simulate .fai file (since samtools faidx is mocked)
    fai_file = output_dir / "all_viruses_dereplicated.fna.fai"
    fai_file.write_text("contig1\t1000\t0\t80\t81\ncontig2\t1000\t0\t80\t81\n")

    # Run
    result = clusterer.cluster(input_contigs, input_reads, output_dir)

    # Verify calls
    assert mock_check_dependency.call_count >= 6

    # Verify commands
    # 1. CD-HIT
    # 2. FastANI
    # 3. Bowtie2 Build
    # 4. CoverM

    # Extract commands passed to subprocess.run
    commands = [call_args[0][0] for call_args in mock_run.call_args_list]

    # Check if key commands are present (by checking executable path mostly which is /usr/bin/tool)
    assert any("cd-hit-est" in str(cmd) or "/usr/bin/tool" in cmd for cmd in commands)
    # Since we mocked check_dependency to return /usr/bin/tool, strict string matching on command name might fail if we used the path
    # But in the code: str(cd_hit) -> "/usr/bin/tool"

    # Let's check descriptions passed to _run_command if we could, but mock_run only sees the list
    # The first arg to run is the command list.

    # We can check if specific flags are present
    assert any("-c" in cmd and "0.99" in cmd for cmd in commands) # CD-HIT
    assert any("--ql" in cmd for cmd in commands) # FastANI
    # Bowtie2 build
    assert any("bowtie2-build" in str(cmd) or "/usr/bin/tool" in cmd for cmd in commands)
    # CoverM
    assert any("coverm" in str(cmd) or "/usr/bin/tool" in cmd for cmd in commands)
    # Samtools faidx
    assert any("faidx" in cmd for cmd in commands)

    # Verify output files exist (mocked creation in _process_clusters)
    assert result.clusters.exists()
    assert result.abundance_table.exists()

    # Verify Abundance Matrix content
    # Since pd.read_csv is mocked, we read the file text directly
    content = result.abundance_table.read_text()
    assert "sample1" in content
    assert "vOTU" in content
