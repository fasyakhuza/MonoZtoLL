My analysis codes to get the dilepton pT distributions after each selection are ```xAna_pre_Zpt_optimization_bkg.C``` and ```xAna_pre_Zpt_optimization_DYincl.C```.

These two macros were submitted to HTCondor and create output root files that I saved in my private eos.

Then, to get the histograms of dilepton pT distributions after each selection, you need to run below commands.

```
git clone https://github.com/fasyakhuza/MonoZtoLL/tree/main/checkDoublePeak.git

mkdir outputplot

root -l xPlot_background_eachType.C

```

The output histograms will be saved in "outputplot" directory
