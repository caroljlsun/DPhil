### ------------------------ ###
# Procedure: Rewrite Peaky filesystem functions
# Project: Capture C
# User: CZJS
# Date: MAR-2025
### ------------------------ ###


library("data.table")
library("tidyverse")
library("rJavaEnv")
library("peaky")
library("purrr")
library("Biobase")

### ------------ Functions to rewrite... ------------ ###


bin_interactions_fs()

model_bin_fs()

split_baits_fs()

peaky_fs()

interpret_peaky_fs()

peaky_wrapup()


### ------------ bin_interactions_fs ------------ ###


bin_interactions_fs_rw <- function (interactions_file, fragments_file, output_dir, min_dist = 2500, 
                                    max_dist = Inf, bins = 5) 
{
    L = paste0(output_dir, "/log_bins.txt")
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    write("BINNING\n", file = L, append = FALSE, sep = "")
    note(L, T, "Reading interactions from ", interactions_file)
    D = fread(interactions_file)
    note(L, T, "Reading fragment information from ", fragments_file)
    fragments = fread(fragments_file)
    BI = bin_interactions(D, fragments, bins, min_dist, max_dist, 
                          L)
    save_bin = function(Dbin, output_dir, bins) {
        filename = paste0(output_dir, "/bin_", as.character(Dbin$dist.bin)[1], 
                          ".rds")
        note(L, F, filename)
        saveRDS(Dbin, file = filename)
    }
    note(L, T, "Saving binned interactions:")
    by(BI$interactions, BI$interactions$dist.bin, FUN = save_bin, 
       output_dir = output_dir, bins = bins)
    note(L, T, "Saving bin details to ", output_dir, "/bins.txt")
    fwrite(BI$bins, file = paste0(output_dir, "/bins.txt"), 
           sep = " ", na = 0)
    return(list(output_dir = output_dir, interactions = BI$interactions, 
                bins = BI$bins))
}





### ------------ model_bin_fs ------------ ###

model_bin_fs_rw <- function (bins_dir, bin_index, output_dir, subsample_size = NA, 
                             gamlss_cycles = 200, gamlss_crit = 0.1, formula_add = NA, 
                             formula_replace = NA) 
{
    L = paste0(output_dir, "/log_fit_", bin_index, ".txt")
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    write("FITTING\n", file = L, append = FALSE, sep = "")
    b = bin_index
    output_path_fit = paste0(output_dir, "/fit_", b, ".rds")
    output_path_residuals = paste0(output_dir, "/residuals_", 
                                   b, ".rds")
    bin_path = paste0(bins_dir, "/bin_", b, ".rds")
    note(L, F, "Loading bin from ", bin_path)
    bin = readRDS(bin_path)
    BM = model_bin(bin, subsample_size = subsample_size, gamlss_cycles = gamlss_cycles, 
                   gamlss_crit = gamlss_crit, formula_add = formula_add, 
                   formula_replace = formula_replace, log_file = L)
    if (is.gamlss(BM$fit)) {
        note(L, F, "Saving fit to ", output_path_fit, "...")
        saveRDS(BM$fit, output_path_fit)
        note(L, T, "Saving all residuals to ", output_path_residuals, 
             "...")
        saveRDS(BM$residuals, output_path_residuals)
    }
    return(list(output_dir = output_dir, fit = BM$fit, residuals = BM$residuals))
}


### ------------split_baits_fs------------ ###


split_baits_fs_rw <- function (bins_dir, residuals_dir, indices, output_dir, plots = TRUE) 
{
    L = paste0(output_dir, "/log_split.txt")
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    write("SPLITTING BAITS\n", file = L, append = FALSE, sep = "")
    bin_files = paste0(bins_dir, "/bin_", indices, ".rds")
    residual_files = paste0(residuals_dir, "/residuals_", indices, 
                            ".rds")
    note(L, T, "Loading bins...")
    note(L, F, bin_files)
    bins = lapply(bin_files, readRDS)
    note(L, T, "Loading residuals...")
    note(L, F, residual_files)
    residuals = lapply(residual_files, readRDS)
    D = split_baits(bins, residuals, log_file = L)
    save_bait = function(baitset, output_dir) {
        filename = paste0(output_dir, "/bait_", as.character(baitset$baitID)[1], 
                          ".rds")
        note(L, F, filename)
        saveRDS(baitset, file = filename)
    }
    plot_bait = function(baitset, output_dir) {
        filename = paste0(output_dir, "/bait_", as.character(baitset$baitID)[1], 
                          ".pdf")
        plot = ggplot(data = baitset, aes(x = dist, y = residual)) + 
            geom_hline(yintercept = 0, color = "grey") + geom_vline(xintercept = 0, 
                                                                    color = "steelblue") + geom_point() + xlab("Distance from bait (bp)") + 
            ylab("Residuals") + ggtitle(paste0("Bait ", as.character(baitset$baitID)[1]))
        ggsave(plot, filename = filename, width = 10, height = 4)
    }
    note(L, T, "Saving baits...")
    by(D, as.factor(D$baitID), FUN = save_bait, output_dir = output_dir)
    bait_paths = paste0(output_dir, "/bait_", unique(D$baitID), 
                        ".rds")
    write(bait_paths, paste0(output_dir, "/baitlist.txt"), sep = "\n")
    if (plots == TRUE) {
        note(L, T, "Plotting baits...")
        by(D, as.factor(D$baitID), FUN = plot_bait, output_dir = output_dir)
    }
    note(L, T, "Done.")
    return(list(output_dir = output_dir, bait_paths = bait_paths))
}


### ------------------------ ###


peaky_fs()
### ------------------------ ###

interpret_peaky_fs()

### ------------------------ ###

peaky_wrapup()



mod_gamlss_model_bin <- function (bin, subsample_size = NA, gamlss_cycles = 200, gamlss_crit = 0.1, 
                                  formula_add = NA, formula_replace = NA, log_file = NA) {
    gamlss.tr::gen.trun(0, family = "DPO", type = "left", name = ".0tr")
    NBI.0tr <<- gamlss.tr::trun(0, family = "NBI", type = "left", 
                                name = ".0tr", local = FALSE)
    L = log_file
    rm(log_file)
    if (is.na(subsample_size)) {
        subsample_size = nrow(bin)
        note(L, T, "No subsampling size provided, using all ", 
             subsample_size, " interactions for the null model regression...")
    }
    else {
        note(L, T, "Subsampling ", subsample_size, "/", nrow(bin), 
             " interactions for the null model regression...")
    }
    subset = sample(1:nrow(bin), subsample_size)
    note(L, T, "Fitting with a maximum of ", gamlss_cycles, 
         " iterations...")
    control = gamlss.control(c.crit = gamlss_crit, n.cyc = gamlss_cycles, mu.step = .5) # modded mu.step
    fit_formula <- "N ~ log(abs(dist)) + b.trans_res  + p.trans_res + sqrt(b.length) + sqrt(p.length)"
    if (!is.na(formula_replace)) {
        fit_formula = formula_replace
    }
    else if (!is.na(formula_add)) {
        fit_formula = paste(fit_formula, formula_add, sep = " + ")
    }
    note(L, F, "Using formula:\n", fit_formula)
    fit = gamlss(as.formula(fit_formula), data = bin[subset, 
    ], family = NBI.0tr, sigma.formula = ~log(abs(dist)), 
    control = control)
    if (is.gamlss(fit)) {
        note(L, T, "Converged: ", ifelse(fit$converged, "YES", 
                                         "NO"), "\nIterations: ", fit$iter, "\n\nCoefficients:")
        note(L, F, paste0(names(coef(fit)), "\t", coef(fit)))
        note(L, T, "Obtaining normalized randomized quantile residuals for the full dataset...")
        residuals_all = gamlss:::rqres(pfun = "pNBI.0tr", type = "Discrete", 
                                       ymin = 1, y = bin$N, mu = predict(fit, what = "mu", 
                                                                         newdata = bin, type = "response", data = bin[subset, 
                                                                         ]), sigma = predict(fit, what = "sigma", newdata = bin, 
                                                                                             type = "response", data = bin[subset, ]))
        return(list(fit = fit, residuals = residuals_all))
    }
    else {
        note(L, T, "No fit object obtained.")
    }
}


models = by(BI_mine$interactions, BI_mine$interactions$dist.bin, mod_gamlss_model_bin, subsample_size=35000)

