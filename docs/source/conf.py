# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('../../python'))

# -- Project information -----------------------------------------------------
project = 'Pysic-rs'
copyright = '2026, HFThot Research Lab'
author = 'HFThot Research Lab'
release = '0.1.0'
version = '0.1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'myst_parser',
    'sphinxcontrib.mermaid',
    'sphinx_copybutton',
    'nbsphinx',
]

# Notebook rendering: notebooks are shipped pre-executed (kernel rhftlab);
# readthedocs must never re-execute them (pysicrs isn't installed there).
nbsphinx_execute = 'never'
nbsphinx_allow_errors = False
nbsphinx_timeout = 300

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', 'notebooks/README.md']

# -- Options for HTML output -------------------------------------------------
html_theme = 'furo'
html_static_path = ['_static']
html_title = 'Pysic-rs Documentation'
html_short_title = 'Pysic-rs'

# Announcement bar (top panel) - Furo supports this via html_theme_options
html_theme_options = {
    "light_css_variables": {
        "color-brand-primary": "#1e3a5f",
        "color-brand-content": "#1e3a5f",
        "color-api-pre-background": "#eef2ff",
        "color-background-primary": "#fefefe",
        "color-background-secondary": "#eef2ff",
        "color-foreground-primary": "#2d2d3a",
        "color-foreground-secondary": "#6b6b7a",
        "color-foreground-muted": "#9a9ab0",
        "color-border": "#c7d2fe",
        "color-announcement-background": "#1e3a5f",
        "color-announcement-text": "#ffffff",
    },
    "dark_css_variables": {
        "color-brand-primary": "#60a5fa",
        "color-brand-content": "#60a5fa",
        "color-api-pre-background": "#0f172a",
        "color-background-primary": "#0f172a",
        "color-background-secondary": "#1e293b",
        "color-foreground-primary": "#e8e8f0",
        "color-foreground-secondary": "#a0a0b8",
        "color-foreground-muted": "#707088",
        "color-border": "#334155",
        "color-announcement-background": "#1e3a5f",
        "color-announcement-text": "#ffffff",
    },
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
    "top_of_page_buttons": ["view", "edit"],
    "source_repository": "https://github.com/ThotDjehuty/pysic-rs",
    "source_branch": "main",
    "source_directory": "docs/source/",
    "announcement": "🚀 <strong>Pysic-rs v0.1.0</strong> — Mathematical physics engine in Rust + Python • textbook-validated • <a href=\"changelog.html\">View Changelog</a> • <a href=\"https://optimiz-rs.readthedocs.io/\">Optimiz-rs Docs</a> • <a href=\"https://hfthot-lab.eu/thotbook-amenti.html\">HFThot Research Lab</a>",
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/ThotDjehuty/pysic-rs",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 24 24">
                    <path fill-rule="evenodd" d="M12 2C6.477 2 2 6.484 2 12.017c0 4.425 2.865 8.18 6.839 9.504.5.092.682-.217.682-.483 0-.237-.008-.868-.013-1.703-2.782.605-3.369-1.343-3.369-1.343-.454-1.158-1.11-1.466-1.11-1.466-.908-.62.069-.608.069-.608 1.003.07 1.531 1.032 1.531 1.032.892 1.53 2.341 1.088 2.91.832.092-.647.35-1.088.636-1.338-2.22-.253-4.555-1.113-4.555-4.951 0-1.093.39-1.988 1.029-2.688-.103-.253-.446-1.272.098-2.65 0 0 .84-.27 2.75 1.026A9.564 9.564 0 0112 6.844c.85.004 1.705.115 2.504.337 1.909-1.296 2.747-1.027 2.747-1.027.546 1.379.202 2.398.1 2.651.64.7 1.028 1.595 1.028 2.688 0 3.848-2.339 4.695-4.566 4.943.359.309.678.92.678 1.855 0 1.338-.012 2.419-.012 2.747 0 .268.18.58.688.482A10.019 10.019 0 0022 12.017C22 6.484 17.522 2 12 2z" clip-rule="evenodd"></path>
                </svg>
            """,
            "class": "",
        },
        {
            "name": "Optimiz-rs",
            "url": "https://optimiz-rs.readthedocs.io/",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 24 24">
                    <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm0 18c-4.41 0-8-3.59-8-8s3.59-8 8-8 8 3.59 8 8-3.59 8-8 8zm-1-13h2v6h-2zm0 8h2v2h-2z"/>
                </svg>
            """,
            "class": "",
        },
        {
            "name": "HFThot Lab",
            "url": "https://hfthot-lab.eu/thotbook-amenti.html",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 24 24">
                    <path d="M12 2L1 21h22L12 2zm0 4.7l6.5 11.3H5.5L12 6.7z"/>
                </svg>
            """,
            "class": "",
        },
    ],
}

# Napoleon settings for Google/NumPy docstring parsing
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = False
napoleon_attr_annotations = True

# Intersphinx configuration
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    # 'optimizrs': ('https://optimiz-rs.readthedocs.io/', None),  # enable once published
}

# MyST parser configuration for markdown support
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "linkify",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]
myst_heading_anchors = 3

# Custom CSS
html_css_files = ["custom.css"]

# Mermaid configuration
mermaid_version = "10.9.0"
mermaid_init_js = "mermaid.initialize({startOnLoad:true, theme:'dark', themeVariables:{primaryColor:'#1e3a5f',primaryTextColor:'#fff',primaryBorderColor:'#1e3a5f',lineColor:'#60a5fa',secondaryColor:'#0f172a',tertiaryColor:'#1e293b'}});"

# Copy button configuration
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True
copybutton_only_copy_prompt_lines = True
copybutton_remove_prompts = True

# Todo extension
todo_include_todos = True
todo_emit_warnings = True
