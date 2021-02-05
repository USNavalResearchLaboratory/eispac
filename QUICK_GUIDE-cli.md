# A Quick and Dirty Guide to the Updated eispac Software Command Line Interface
_version 0.7_

These are helper routines to get you started quickly.

## eis_download_files

Download data files from https://eis.nrl.navy.mil

```
> eis_download_files eis_l0_20201227_095348
```

## eis_browse_templates

A widget that allows you to select a file and see what fitting templates are available. First click
on "Select a header file," navigate to the file you just downloaded, select it. A good place to
start is with a single component fit of Fe XII 195.119. Select that from the list and click "Copy
Template to Local Dir". 

```
> eis_browse_templates
```

After running this you should have a new subdirectory `eis_template_dir` with the file
`fe_12_195_119.1c.template.h5` in it.

## eis_fit_files

Given some data files and some templates perform some fits. Inputs can be file names or directories.

```
> eis_fit_files data_eis/ eis_template_dir/
```

## eis_plot_fit

```
> eis_plot_fit data_eis/eis_20201227_095348.fe_12_195_119.1c-0.fit.h5
```

