""" Development version of the EIS Python Analysis Code """

import setuptools

setuptools.setup(
    name = 'eispac',
    version = '0.8.0',
    description = 'Python analysis tools for Hinode / EIS data',
    author = 'NRL EISPAC Development Team',
    author_email = 'N/A',
    license = 'MIT',
    url = "https://github.com/hpwarren/pyEIS-test",
    project_urls = {
        "Data": "https://eis.nrl.navy.mil/"},
    packages = setuptools.find_packages(),
    entry_points = {
        "console_scripts": [
            "eis_browse_templates = scripts.eis_browse_templates:eis_browse_templates",
            "eis_catalog = scripts.eis_catalog:eis_catalog",
            "eis_download_files = scripts.eis_download_files:eis_download_files",
            "eis_fit_files = scripts.eis_fit_files:eis_fit_files",
            "eis_plot_fit = scripts.eis_plot_fit:eis_plot_fit"]
        },
    classifiers = [
        "Development Status :: 3 - Alpha"
        "Programming Language :: Python :: 3.7",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"],
    keywords = 'solar, sun, physics, spectroscopy, Hinode, EIS',
    python_requires = '>=3.7',
    install_requires = [
        "numpy",
        "scipy",
        "matplotlib",
        "h5py",
        "astropy",
        "sunpy",
        "ndcube",
        "wget"],
    include_package_data = True
    )
