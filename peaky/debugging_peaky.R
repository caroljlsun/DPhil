bin_interactions_mine = function (interactions, fragments, bins = 5, min_dist = 2500, 
          max_dist = Inf, log_file = NA) 
{
    
    interactions = interactions_mine
    fragments = fragments_mine
    bins = 5
    min_dist = 2500
    max_dist = Inf
    log_file = NA
    
    
    
    D = interactions
    rm(interactions)
    L = log_file
    rm(log_file)
    if (ncol(D) < 3 | all.equal(colnames(D)[1:3], c("baitID", 
                                                    "preyID", "N")) != TRUE) {
        note(L, T, "First three columns of interactions file should be named: baitID preyID N")
        stop()
    }
    if (ncol(fragments) < 4 | !(all(colnames(fragments)[1:4] == 
                                    c("chrom", "chromStart", "chromEnd", "ID")) | all(colnames(fragments)[1:4] == 
                                                                                      c("chr", "start", "end", "ID")))) {
        note(L, T, "First four columns of fragments file should be named: chrom chromStart chromEnd ID")
        stop()
    }
    data.table::setnames(fragments, c("chr", "start", "end", "ID"))
    note(L, T, "Calculating fragment characteristics...")
    fragments = fragments[, .(chr, mid = round((end + start)/2), 
                              length = end - start, ID)]
    note(L, T, "Adding fragment characteristics for baits...")
    D = merge(D, fragments[fragments$ID %in% D$baitID], by.x = "baitID", 
              by.y = "ID", all.x = TRUE)
    data.table::setnames(D, paste0(ifelse(names(D) %in% names(fragments), 
                              "b.", ""), names(D)))
    note(L, T, "Adding fragment characteristics for preys...")
    D = merge(D, fragments[fragments$ID %in% D$preyID], by.x = "preyID", 
              by.y = "ID", all.x = TRUE)
    data.table::setnames(D, paste0(ifelse(names(D) %in% names(fragments), 
                              "p.", ""), names(D)))
    note(L, T, "Calculating interaction distances...")
    D[b.chr == p.chr, `:=`(dist, (p.mid - b.mid))]
    note(L, T, "Calculating total trans-chromosomal read counts for each bait...")
    trans = D[, .(b.trans = sum((b.chr != p.chr) * N)), by = .(baitID, 
                                                               b.chr)]
    note(L, T, "Modelling those as a function of bait chromosome...")
    if (length(unique(D$b.chr)) == 1) {
        note(L, T, "All baits actually appear to be on the same chromosome, reverting to an intercept only model...")
        trans_model = gamlss(b.trans ~ 1, data = trans, family = gamlss.dist::NBI)
    } else {
        trans_model = gamlss(b.trans ~ as.factor(b.chr), data = trans,
                             family = gamlss.dist::NBI)
        # trans_model = gamlss(b.trans ~ as.factor(b.chr), data = na.omit(trans), family = "DPO")
    }
    trans[, `:=`(b.trans_res, trans_model$residuals)]
    head(trans)
    D = merge(D, trans[, .(baitID, b.trans, b.trans_res)], by = "baitID", 
              all.x = TRUE)
    note(L, T, "Adding trans-chromosomal interactivity covariate for preys that were also baited (0 for preys not baited)...")
    D[, `:=`(p.trans_res, b.trans_res[match(preyID, baitID)])][is.na(p.trans_res), 
                                                               `:=`(p.trans_res, 0)]
    note(L, T, "Excluding ", D[, sum(abs(dist) < min_dist, na.rm = TRUE)], 
         " interactions that are too proximal (distance < ", 
         min_dist, " bp)...")
    D = D[is.na(dist) | abs(dist) >= min_dist]
    note(L, T, "Excluding ", D[, sum(abs(dist) > max_dist, na.rm = TRUE)], 
         " interactions that are too distal (distance > ", max_dist, 
         " bp)...")
    D = D[is.na(dist) | abs(dist) <= max_dist]
    note(L, T, "Assigning ", bins, " distance bins...")
    D[p.chr == b.chr, `:=`(dist.bin, as.factor(as.integer(cut_number(D[p.chr == 
                                                                           b.chr, abs(dist)], n = bins))))]
    overview = D[, .(dist.abs.min = min(abs(dist)), dist.abs.max = max(abs(dist)), 
                     interactions = .N), by = dist.bin][order(dist.abs.min), 
                     ]
    note(L, T, "Done.")
    return(list(interactions = D, bins = overview))
}
