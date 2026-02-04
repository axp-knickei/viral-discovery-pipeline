import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.abspath("../../src"))

project = "Viral Discovery Pipeline"
author = "Viral Pipeline Team"
release = "0.1.0"

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
}

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
]

templates_path = ["_templates"]
exclude_patterns = []

html_theme = "alabaster"
html_last_updated_fmt = "%Y-%m-%d"
html_theme_options = {
    "description": f"Last updated {datetime.utcnow().date()}",
}
