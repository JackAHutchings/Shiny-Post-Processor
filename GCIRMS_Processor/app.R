# Version 0!

rm(list=ls())

library(dplyr) # Data processing
library(ggplot2); theme_set(theme_classic()) # Data visualization
library(tidyr) # Data processing
library(zoo) # needed for na.locf function
library(lemon) # needed for facet_rep_wrap/grid functions
library(writexl) # export the finished data
library(readxl) # import the template file
library(shinydashboard) # shiny
library(httr) # Needed for fetching URLs
library(DT) # tabular data is rendered via DT
library(shiny) # shiny

#testing
# {
#     raw_isodat_file = "(H2 Export).csv"
#     gcirms_template_file = "GCIRMS Template - Phthalic Acid.xlsx"
#     compound_option = "Assigned by IRMS export in Component/comp column."
#     drift_option = "No drift correction (use this when there is no apparent drift or if you use bracketed scale normalization)"
#     drift_comp = "Mean drift of all compounds"
#     size_option = "Peak Height (amplitude, mV)"
#     size_cutoff = 3300
#     size_normal_peak_action = "Remove size effect with 'Normal' size effect function."
#     size_small_peak_action = "Remove size effect with 'Small' size effect function."
#     size_toosmall_peak_action = "Remove these from results."
#     size_large_peak_action = "No size effect correction."
#     acceptable_peak_units = "Peak Height (amplitude, mV)"
#     largest_acceptable_peak = NA
#     smallest_acceptable_peak = NA
#     normalization_option = "Linear interpolation between adjacent normalization standards"
#     normalization_comps = c("C20 FAME", "C31 Alkane")
#     derivatization_option = "Template-defined derivative \u03B4D."
# 
# }

ingest_function <- function(raw_isodat_file,gcirms_template_file) {
    
    # Read in the 'raw' IRMS export}
    if(grepl("csv",raw_isodat_file)){raw_irms_data <- read.csv(raw_isodat_file)} else 
    if(grepl("xlsx",raw_isodat_file)){raw_irms_data <- read_excel(raw_isodat_file)}
    
    sample_info <- read_xlsx(gcirms_template_file,sheet="Samples",range=cell_cols("A:G")) %>%  # Read in the 'Samples' sheet
        rename(id1 = `Identifier 1`) # Change to the variable used throughout this script...
    standard_info <- read_xlsx(gcirms_template_file,sheet="Standards",range=cell_cols("A:G")) %>% rename(std_mix = `Identifier 1`) # Read in the 'Standards' sheet
    headers <- read_xlsx(gcirms_template_file,sheet="Headers",range=cell_cols("B"))$export_header %>% # Vector of the headers in the raw_irms_data
        setNames(read_xlsx(gcirms_template_file,sheet="Headers",range=cell_cols("A"))$ingest_header) # Values to rename the original headers
    derivatization_table <- read_xlsx(gcirms_template_file,sheet="Derivatization",range=cell_cols("A:F")) # Read in the 'Derivatization' Sheet
    initials <- read_xlsx(gcirms_template_file,sheet="Initials",range=cell_cols(c("A:B"))) %>% # Read in the 'Initials' sheet
        group_by(variable) %>% mutate(observation = 1:n()) %>% spread(variable,initial_value) # Reformat so that each initial is a variable.
    rt_table <- read_xlsx(gcirms_template_file,sheet="Retention Times",range=cell_cols("A:C")) %>% 
        rename(rt_table = rt,comp_rt_table = comp_class)
    
    header_mismatch <- names(headers[which(!(headers %in% colnames(raw_irms_data)))]) # A list of initially missing headers.
    optional_headers <- c("raw_R","conc","peak","rt") # A list of optional headers.
    missing_headers <- header_mismatch[which(!(header_mismatch %in% optional_headers))] # A list of non-optional missing headers.
    headers_to_use <- headers[which(!(names(headers) %in% header_mismatch))] # Original header list sans *any* missing headers.
    
    # If there are no non-optional missing headers, then go ahead with initial data processing:
    if(length(missing_headers) == 0){
        #Data Ingestion. Yum!
        {
        rawdata <- raw_irms_data %>% # Reads in the file.
            filter(!is.na(Row)) %>% # Remove any blank rows (sometimes caused by export weirdness).
            rename(!!!headers_to_use) %>%  # Rename the headers we will use based on the 'Headers' sheet in the GCIRMS template.
            mutate(row = row - min(row) + 1) %>%  # Generates a row index starting with first exported injection. This is helpful as the first few rows are usually unused warmups.
            rowwise() %>% 
            mutate(rt = ifelse("rt" %in% names(.),rt,NA)) %>% # Check if retention time was provided. If not, enter NA. RT is only used to assign compound names.
            mutate(comp = ifelse("comp" %in% names(.),as.character(comp),NA)) %>% # Check if compound names are provided. If not, enter NA as a placeholder.
            mutate(peak = ifelse("peak" %in% names(.),peak,1:n())) %>% # Check if the peak column is present. If it isn't, generate an index so the code still functions.
            mutate(conc = ifelse("conc" %in% names(.),conc,1)) %>% # Check if concentrations were provided. If not, use '1' as the placeholder. In principle, the
            # end-user could provide unique Identifier 1 entries for each concentration of a standard and have those match the IRMS data. This should be done if
            # the end-user has size effect standards with varying compound concentration ratios (i.e., if each level were manually mixed rather than serially diluted)
            mutate(raw_R = ifelse("raw_R" %in% names(.),raw_R,0)) %>% #Check if the raw element ratio was provided. If not, enter zero. The raw ratio is only used as
            # a diagnostic plot.
            as.data.frame() %>% # This 'turns off' the rowwise mode.
            group_by(row,id1,id2,comp) %>% 
            filter(peak == max(peak)) %>% # Manually added peaks are added to the end of the peak list, but with the same 'comp' value. This filters
            # so that only the last peak of a comp is used. This behavior is only known to be true for Isodat 3.0! Older/other software may have distinct behavior.
            # The safest way to ensure compatibility is to just have a single observation for each compound in each injection.
            filter(comp != "" & comp != "-") %>% # Remove non-identified peaks.
            select(row,id1,id2,conc,comp,rt,area,ampl,raw_R,dD_raw) %>% # Select only the headers we will use.
            arrange(row) %>% ungroup()
        
        injections <- rawdata %>% 
            select(row,id1) %>% #row is our unique injection identifier, but we do multiple injections per vial!
            distinct() %>% # Removes all the extra rows for the exported compound information.
            arrange(row) %>% # Sort the data frame by row.
            group_by(id1) %>% # For each unique id1 we'll do the following... (warning: janky solution)
            mutate(row_diff = c(0,diff(row)), # What's the difference in row for each injection of a vial? Practically, all non-first injections in each series is equal to 1.
                   row_base = ifelse(row_diff == 1,NA,row), # If an injection isn't the first in a replicate injection series, mark it NA.
                   row_base = na.locf(row_base)) %>% # locf = last observation carried forward, this marks each subsequent injection in a series with its starting row.
            select(-c(row_diff))
        
        compounds <- rawdata %>% mutate(dummy=T) %>% 
            full_join(rt_table %>% mutate(dummy=T),by="dummy") %>% 
            group_by(row,rt) %>% 
            mutate(rt_diff = abs(rt - rt_table)) %>% 
            group_by(rt_table,window) %>% 
            filter(rt_diff < window) %>% 
            group_by(row,comp_rt_table) %>% 
            mutate(closest_RT = ifelse(rt_diff == min(rt_diff),T,F)) %>% 
            filter(closest_RT) %>% 
            select(row,id1,id2,rt,comp_rt_table)
        
        output <- full_join(rawdata,injections,by = c("row", "id1")) %>% 
            rename(comp_irms_output = comp) %>% 
            full_join(compounds,by=c("row","id1","id2","rt"))
        }
    
        list("raw_irms_check" = T,
             "raw_irms_data" = raw_irms_data,
             "sample_info" = sample_info,
             "standard_info" = standard_info,
             "derivatization_table" = derivatization_table,
             "initials" = initials,
             "rt_table" = rt_table,
             "output" = output)
    } else # Otherwise, check for missing headers and list them as the only output.
    if(length(missing_headers) > 0){list("raw_irms_check" = missing_headers)}
}

comp_assign_function <- function(input,standard_info,compound_option) {
    
    compound_action <- ifelse(compound_option == "Assigned by IRMS export in Component/comp column.","comp_irms_output",
                       ifelse(compound_option == "Assigned by retention time matching with Retention Times sheet.","comp_rt_table",
                              "compound_option in comp_assign_function broken! Check that the radioButton choices still match!"))
    
    data <- input %>% 
        gather(compound_option,comp,comp_irms_output,comp_rt_table) %>% 
        filter(compound_option %in% compound_action) %>% 
        distinct() %>% 
        select(-compound_option) %>% 
        mutate(id1 = as.character(id1),
               comp = as.character(comp)) %>% 
        mutate(std_mix = ifelse(id2 == "sample" | id2 == "derivatization",NA,id1)) %>% #id1 should be used to identify the standard mix names.
        separate(comp,c("comp","class"),sep=" ",fill="right") %>% 
        left_join(standard_info,by = c("std_mix", "comp", "class")) %>%  # Join with standard_info
        mutate(conc = conc * relative_concentration) %>% 
        filter(!is.na(row)) %>% 
        mutate(dD_processing = dD_raw) %>%  # We use the 'dD_processing' variable as the placeholder for each data correction step.
        unite(comp_class,c(comp,class),sep=" ",remove=F,na.rm=T) # This code currently assumes that compound and class uniquely identifies standards in a mix.

    #IRMS Drift Plot
    irms_drift_calc <- data %>%
        mutate(slope = lm(dD_raw~raw_R)$coefficients[2], # Figure out the general peak area ratio to dD relationship... slope first...
               intercept = lm(dD_raw~raw_R)$coefficients[1], # ... then intercept
               rawraw_dD = raw_R * slope + intercept) %>%  # Estimate the dD of each ref peak using the above relationship.
        filter(grepl("Ref",comp)) %>%
        mutate(leftcenter_dD = rawraw_dD - mean(rawraw_dD[which(row==min(row))])) %>% # Center the extra raw 'rawraw_dD' on the first injection, so that we can see how drift proceeded.
        group_by(row) %>%
        mutate(sd = sd(leftcenter_dD)) %>% ungroup() %>% mutate(mean_sd = round(mean(sd),3), sd_sd = round(sd(sd),3)) # Summary stats on how the ref peaks performed.
    
    if(length(irms_drift_calc$comp) == 0){
        plot_irms_drift <- ggplot()+annotate("text",x=0,y=0,label="No Reference Gas Peaks Found.\n\nReference gas peaks should be identified and their name should start with\n'Ref' (i.e., Ref 1 or Reference 1) in the Component/comp column.",color="red")+theme_void()
    } else {
        # Note! Because each injection is calculated relative to its reference peak, this drift is fully incorporated into the "dD_raw" value that we calibrate.
        # As such, we do not need to correct for this drift at all. Instead, it is just a helpful visualization of how much drift the IRMS experienced during the run.
        # If this is huge (i.e., > 25 permil in range), then maybe it can be diagnostic of a problem...
        plot_irms_drift <- ggplot(irms_drift_calc,aes(x=row,y=leftcenter_dD)) +
            geom_point() +
            labs(title = paste0("Mean Per-Injection Ref Peak SD: ",irms_drift_calc$mean_sd[1]," +/- ",irms_drift_calc$sd_sd[1],"\u2030"),
                 x = "Injection # in Sequence",
                 y = "IRMS Drift (\u2030)")    
    }

    #Area/Amplitude Curves
    area_ampl_curve_data <- data %>% filter(grepl("size",id2) & comp != "Ref") # Grab just the size effect standards
    
    if(length(area_ampl_curve_data$id1) == 0) {
        plot_area_curve <- ggplot()+annotate("text",x=0,y=0,label="No Size Effect Standards Found.",color="red")+theme_void()
        plot_ampl_curve <- ggplot()+annotate("text",x=0,y=0,label="No Size Effect Standards Found.",color="red")+theme_void()
    } else {
        plot_area_curve <- ggplot(area_ampl_curve_data,aes(x=conc,y=area)) +
            geom_point() +
            geom_smooth(method="lm",se=F) +
            facet_rep_wrap(~class+comp,nrow=2) +
            labs(x = "Concentration (ng / uL)",
                 y = "Peak Area")
        plot_ampl_curve <- ggplot(area_ampl_curve_data,aes(x=conc,y=ampl)) +
            geom_point() +
            geom_smooth(method="lm",se=F) +
            facet_rep_wrap(~class+comp,nrow=2) +
            labs(x = "Concentration (ng / uL)",
                 y = "Peak Amplitude")
    }
        
    list("data" = data,                 
         "plot_irms_drift" = plot_irms_drift,
         "plot_area_curve" = plot_area_curve,
         "plot_ampl_curve" = plot_ampl_curve)
}

drift_function <- function(input,drift_option,drift_comp) {
    
    drift_raw <- input %>% filter(comp != "Ref") %>%  # Removes the reference peaks.
        filter(grepl("drift",id2)) %>%   # drift identification may be in id2 or id3
        filter(!is.na(dD_known))
    
    if(length(drift_raw$id2) > 0) {
        drift_raw_plot <-   ggplot(drift_raw,aes(x=row,y=dD_processing)) +
            geom_point() +
            geom_smooth(method="lm",se=F) +
            facet_wrap(~comp_class,scales="free_y") +
            labs(title = "Uncorrected Drift Standards",
                 x = "Injection # in Sequence",
                 y = "Raw \u03B4D (\u2030)")
        drift_comp <- ifelse(drift_comp == "Mean drift of all compounds",NA,drift_comp)
        drift_type <- ifelse(drift_option == "Linear interpolation between adjacent drift samples (use this when drift appears non-linear)",1,
                      ifelse(drift_option == "Linear regression across all drift samples (use this when drift appears linear)",2,
                      ifelse(drift_option == "No drift correction (use this when there is no apparent drift or if you use bracketed scale normalization)",3,
                             "drift_option in drift_function broken! Check that the radioButton choices still match!")))
        
        drift_correction <- drift_raw %>%
            group_by(comp_class) %>% 
            mutate(i = dense_rank((row_base)),
                   dD_processing = dD_processing - mean(dD_processing)) %>%  # Generate an indexing variable for each set of standard injections... only relevant using drift_type 1
            filter(ifelse(is.na(drift_comp),T,comp_class==drift_comp)) %>% 
            group_by(row,i,row_base) %>% 
            summarize(dD_mean = mean(dD_processing,na.rm=T)) %>% 
            group_by(i) %>%
            # The next mutate is doing the heavy lifting of calculating the slopes. Each 'bin' is a given drift sample (i) and the next drift sample (i+1)
            # Since injection-to-injection time is practically constant, row number (row_mean) takes the place of time when calculating drift.
            # We are basically just doing a two point linear regression to find the slope and y-intercept of the line.
            mutate(slope_1 = lm(dD_mean~row,data=na.omit(.[which((.$i == i[1] | .$i == i[1]+1 )),]))$coefficients[2],
                   intercept_1 = lm(dD_mean~row,data=na.omit(.[which((.$i == i[1] | .$i == i[1]+1 )),]))$coefficients[1]) %>% 
            ungroup() %>% 
            mutate(slope_1 = ifelse(i == max(i),slope_1[which(i == max(i)-1)],slope_1), # Assigns second to last slope to the last drift set.
                   intercept_1 = ifelse(i == max(i),intercept_1[which(i == max(i)-1)],intercept_1)) %>% # Same as above.
            # Next, we need to make a data frame to fill in all the rows (injections) that occurred between the drift samples.
            full_join(data.frame(row = seq(1,max(input$row))),by="row") %>% 
            arrange(row) %>%  # Sort the dataset by compound and then by row. This is important because we use a 'last-observation-carried-forward' command...
            mutate(slope_2 = lm(dD_mean~row)$coefficients[2], # Calculates the slope for drift_type = 2
                   intercept_2 = lm(dD_mean~row)$coefficients[1], # Calculates the intercept for drift_type = 2
                   slope_1 = na.locf(slope_1), # The slope_1 and slope_1 have been calculated for each drift standard, but we want to use that standard's
                   intercept_1 = na.locf(intercept_1)) %>% # slope and intercept for every injection until the next drift standard. The na.locf command accomplishes this.
            # Okay, so correction_1 needs some exposition. Because we are doing linear interpolation between each nearest pair of drift standards, we are predicting
            # the dD of a drift standard between those two within the 'round()' command of correction_1. Then, we subtract the grand mean (dD_center) from that value
            # to determine the absolute dD deviation of a drift standard at that row. This estimates the absolute distance (in permil) that a drift standard in a given
            # row would be from the grand mean of dD of actually injected drift standards. The assumption is, of course, that the drift between adjacent drift standards
            # is linear.
            mutate(correction_1 = round(row * slope_1 + intercept_1,6),
                   correction_2 = (row * slope_2 + intercept_2), # This is just linear regression across all the drift samples.
                   correction_3 = 0) %>% # No correction.
            gather(correction_method,drift_correction,correction_1,correction_2,correction_3) %>% #Collapse all of the possible corrections into two columns...
            separate(correction_method,c("var","method")) %>% # Separate the correction identifier...
            filter(method %in% drift_type) #... and filter the data set for the correction indicated in the drift_type object.
        
        # This applies the drift correction to the data set.
        output <- full_join(input %>% filter(comp != "Ref"), # Removes the reference peaks.
                                      drift_correction %>% select(row,method,drift_correction) %>% distinct(), #Cleans up the drift_correction a bit.
                                      by = "row") %>% 
            mutate(dD_predrift = dD_processing,
                   dD_processing = dD_processing - drift_correction,
                   dD_drift = dD_processing)
        
        drift_corrected_standards_plot <- ggplot(output %>% filter(grepl("drift",id2)),aes(x=row,y=dD_processing)) +
            geom_segment(aes(x = row, xend = row, y = dD_predrift, yend = dD_processing)) +
            geom_point(color="blue",size=1.5) +
            geom_point(aes(x=row,y=dD_predrift),size=1.5) +
            geom_smooth(method="lm",se=F) +
            facet_wrap(~comp,scales="free_y") +
            labs(title = "Drift-Corrected Drift Standards",
                 subtitle = "Blue = Corrected, Black = Uncorrected",
                 x = "Injection # in Sequence",
                 y = "Uncalibrated \u03B4D (\u2030)")
        
        drift_correction_by_row_plot <- ggplot(drift_correction,aes(x=row,y=drift_correction)) +
            geom_hline(aes(yintercept = 0),color="blue") +
            geom_segment(aes(x = row, xend = row, y = drift_correction, yend = 0),arrow = arrow(length=unit(2,"mm"))) +
            geom_point() +
            labs(title = "Absolute drift correction applied to each row.",
                 x = "Injection # in Seqence",
                 y = "Drift Correction (\u2030)")
        
        list("output" = output,
             "drift_correction" = drift_correction,
             "drift_raw_plot" = drift_raw_plot,
             "drift_corrected_standards_plot" = drift_corrected_standards_plot,
             "drift_correction_by_row_plot" = drift_correction_by_row_plot)
    } else {
        output <- input %>% filter(comp != "Ref") %>% 
            mutate(dD_predrift = dD_processing,
                   dD_processing = dD_predrift,
                   dD_drift = dD_processing)
        
        list("output" = output,
             "drift_correction" = ggplot()+annotate("text",x=0,y=0,label="No Drift Standards Found.",color="red")+theme_void(),
             "drift_raw_plot" = ggplot()+annotate("text",x=0,y=0,label="No Drift Standards Found.",color="red")+theme_void(),
             "drift_corrected_standards_plot" = ggplot()+annotate("text",x=0,y=0,label="No Drift Standards Found.",color="red")+theme_void(),
             "drift_correction_by_row_plot" = ggplot()+annotate("text",x=0,y=0,label="No Drift Standards Found.",color="red")+theme_void())
    }
    
}

size_function <- function(input,size_option,size_cutoff,size_normal_peak_action,size_small_peak_action,
                          size_toosmall_peak_action,size_large_peak_action,acceptable_peak_units,
                          largest_acceptable_peak,smallest_acceptable_peak) {

    size_filter <- ifelse(size_option == "Peak Height (amplitude, mV)","ampl",
                   ifelse(size_option == "Peak Area (Vs)","area",
                   ifelse(size_option == "No size effect correction.","ampl",
                    "size_option in size_function broken! Check that the radioButton choices still match!")))
    
    size_opt_out <- ifelse(size_option == "No size effect correction.",0,1)
    
    size_label = factor(size_filter,levels = c("area","ampl"),labels = c("Area ( Vs )","Amplitude ( mV )"))
    
    size_normal_peak_action = ifelse(size_normal_peak_action == "No size effect correction.",0,1)
    size_small_peak_action = ifelse(size_small_peak_action == "No size effect correction.",0,1)
    size_large_peak_action = ifelse(size_large_peak_action == "Remove size effect with 'Normal' size effect function.",1,
                             ifelse(size_large_peak_action == "No size effect correction.",0,
                             ifelse(size_large_peak_action == "Remove these from results.",NA,
                                    "size_large_peak_action in size_function broken! Check that the radioButton choices still match!")))
    size_toosmall_peak_action = ifelse(size_toosmall_peak_action == "Use 'Small' size effect function.",1,
                                ifelse(size_toosmall_peak_action == "No size effect correction.",0,
                                ifelse(size_toosmall_peak_action == "Remove these from results.",NA,
                                       "size_toosmall_peak_action in size_function broken! Check that the radioButton choices still match!")))
    
    acceptable_peak_units = ifelse(acceptable_peak_units == "Peak Height (amplitude, mV)","ampl",
                            ifelse(acceptable_peak_units == "Peak Area (Vs)","area",
                                   "acceptable_peak_units in size_function broken! Check that the radioButton choices still match!"))
    smallest_acceptable_peak = ifelse(is.na(smallest_acceptable_peak),0,smallest_acceptable_peak)
    largest_acceptable_peak = ifelse(is.na(largest_acceptable_peak),Inf,largest_acceptable_peak)
    
    size_effect_input <- input %>% filter(grepl("size",id2))
    
        if(length(size_effect_input$id2) > 0) {
            size_effect_raw <- input %>% filter(comp != "Ref") %>% 
                filter(grepl("size",id2)) %>% # Grab just the size effect standards
                gather(size_var,value,area,ampl) %>% group_by(row,id1,comp,class) %>% 
                mutate(acceptable_peak_value = value[which(size_var == acceptable_peak_units)]) %>% 
                filter(acceptable_peak_value > smallest_acceptable_peak & acceptable_peak_value < largest_acceptable_peak) %>% 
                filter(size_var == size_filter) %>% 
                mutate(size_group = ifelse(value < size_cutoff,"Small","Normal")) %>% 
                group_by(comp) %>% 
                mutate(dD_zeroed = dD_processing - min(dD_processing)) %>% # Zero-centered. We are grouped by comp here so that each compound is zero centered by its mean.
                group_by(size_group) %>% 
                mutate(size_slope = lm(dD_zeroed~value)$coefficients[2],  # Calculate the size effect using the 'value' variable
                       size_intercept = ifelse(!is.na(size_group),lm(dD_zeroed~value)$coefficients[1],0)) %>%  # Doing segmented size correction requires an intercept.
                mutate(size_slope = size_slope * size_opt_out, # Change the size effect slope to zero if we aren't using it.
                       size_intercept = size_intercept * size_opt_out) %>%  # Same for intercept.
                mutate(size_slope = ifelse(grepl("Normal",size_group),size_slope*size_normal_peak_action,size_slope),
                       size_intercept = ifelse(grepl("Normal",size_group),size_intercept*size_normal_peak_action,size_intercept)) %>% 
                mutate(size_slope = ifelse(grepl("Small",size_group),size_slope*size_small_peak_action,size_slope),
                       size_intercept = ifelse(grepl("Small",size_group),size_intercept*size_small_peak_action,size_intercept)) %>% 
                mutate(dD_processing = dD_processing - (size_slope*value + size_intercept)) %>%  # Perform the correction by subtracting the dD 'anomaly'.)
                group_by(comp) %>% 
                mutate(dD_drift_size_zeroed = dD_processing - mean(dD_processing)) %>% 
                ungroup() %>% 
                mutate(size_upper = max(value),
                       size_lower = min(value),
                       size_group = ifelse(is.na(size_group),"Normal",size_group))
            
            size_raw_bycomp_plot <- size_effect_raw %>% ggplot(aes(x=value,y=dD_zeroed)) +
                geom_point() +
                facet_wrap(~comp,scales="free",nrow=1) +
                geom_smooth(method="lm",se=F) +
                labs(title="Uncorrected plot of size effect standards.\nFacetted by compound.",
                     x=size_label,
                     y="Uncorrected \u03B4D ( \u2030 )")
            
            size_raw_grouped_plot <- size_effect_raw %>% ggplot(aes(x=value,y=dD_zeroed,color=size_group)) +
                geom_point() +
                geom_smooth(method="lm",se=F) +
                labs(title="Uncorrected plot of size effect standards.\nThe slope of this line is used for correction.",
                     y="Uncorrected \u03B4D ( \u2030 )",
                     x=size_label)
            
            size_corrected_bycomp_plot <- size_effect_raw %>% ggplot(aes(x=value,y=dD_processing)) +
                geom_point() +
                facet_wrap(~comp,scales="free_y",nrow=1) +
                geom_smooth(method="lm",se=F) +
                labs(title="Corrected plot of size effect standards.\nNote: Slope is not zero because the compounds do not behave 100% identically.",
                     y="Size Corrected \u03B4D ( \u2030 )",
                     x=size_label)
            
            size_corrected_grouped_plot <- size_effect_raw %>% ggplot(aes(x=value,y=dD_drift_size_zeroed)) +
                geom_point(aes(color=comp)) +
                geom_smooth(method="lm",se=F) +
                labs(title="Corrected plot of size effect standards.",
                     y="Size Corrected \u03B4D ( \u2030 )",
                     x=size_label)
            
            size_effect_table <- data.frame(size_group = c("Too Small","Small","Normal","Large")) %>% mutate(size_group = as.character(size_group)) %>% 
                left_join(size_effect_raw %>% select(size_group,size_slope,size_intercept) %>% distinct(),by="size_group") %>% 
                gather(size_coef,value,size_slope,size_intercept) %>% 
                spread(size_group,value) %>% group_by(size_coef) %>% 
                mutate(Small = ifelse(is.na(Small),Normal,Small),
                       `Too Small` = Small*size_toosmall_peak_action,
                       Large = Normal*size_large_peak_action) %>% 
                gather(size_group,value,Large,Normal,Small,`Too Small`) %>% 
                spread(size_coef,value)
            
            output <- input %>% filter(comp != "Ref") %>% 
                gather(size_var,value,area,ampl) %>% 
                group_by(row,comp,class) %>% 
                mutate(acceptable_peak_value = value[which(size_var == acceptable_peak_units)]) %>% 
                filter(acceptable_peak_value > smallest_acceptable_peak & acceptable_peak_value < largest_acceptable_peak) %>% 
                mutate(size_value = value[which(size_var == size_filter)],
                       size_lower = size_effect_raw$size_lower[1],
                       size_upper = size_effect_raw$size_upper[1],
                       size_group = ifelse(size_value < size_cutoff & size_value >=size_lower,"Small","Normal"),
                       size_group = ifelse(is.na(size_group),"Normal",size_group),
                       size_group = ifelse(size_value < size_lower, "Too Small",
                                    ifelse(size_value > size_upper, "Large",size_group))) %>% 
                spread(size_var,value) %>% 
                left_join(size_effect_table,by="size_group") %>% # Get size effect coefficients...
                mutate(dD_presize = dD_processing,
                       dD_processing = dD_presize - (size_slope*size_value + size_intercept),
                       dD_size = dD_processing) %>% 
                filter(!is.na(dD_processing))
            
            size_correction_plot <- data.frame(variable = size_filter,
                                                 size_value = seq(from=0,to=max(output$size_value),by=100),
                                                 size_lower = size_effect_raw$size_lower[1],
                                                 size_upper = size_effect_raw$size_upper[1]) %>% 
                mutate(size_group = ifelse(size_value < size_cutoff & size_value >=size_lower,"Small","Normal"),
                       size_group = ifelse(is.na(size_group),"Normal",size_group),
                       size_group = ifelse(size_value < size_lower, "Too Small",
                                           ifelse(size_value > size_upper, "Large",size_group))) %>% 
                full_join(size_effect_table,by="size_group") %>% 
                mutate(size_correction = (size_value*size_slope+size_intercept),
                       size_correction = ifelse(is.na(size_correction),0,size_correction),
                       corrected = 0) %>% 
                mutate(size_value = ifelse(size_group == "Too Small" & is.na(size_toosmall_peak_action),NA,size_value)) %>% 
                filter(!is.na(size_value)) %>% 
                mutate(size_group = factor(size_group,levels=c("Too Small","Small","Normal","Large"),ordered=T,
                                           labels=c("Smaller Than Size Effect Curve",
                                                    "Lower Segment of Size Effect Curve",
                                                    "Upper Segment/Whole of Size Effect Curve",
                                                    "Larger Than Size Effect Curve"))) %>% 
                ggplot(aes(x=size_value,y=size_correction,color=size_group)) +
                geom_point() +
                geom_segment(aes(x = size_value, xend = size_value, y = size_correction, yend = corrected,group = size_value)) +
                scale_color_brewer(type="qual",palette=2) +
                labs(x=size_label,
                     y="Size Effect Correction (\u2030)",
                     color="Size Group")
                
                
            list("output" = output,
                 "size_filter" = size_filter,
                 "size_effect_raw" = size_effect_raw,
                 "size_raw_bycomp_plot" = size_raw_bycomp_plot,
                 "size_raw_grouped_plot" = size_raw_grouped_plot,
                 "size_corrected_bycomp_plot" = size_corrected_bycomp_plot,
                 "size_corrected_grouped_plot" = size_corrected_grouped_plot,
                 "size_correction_plot" = size_correction_plot
                 )
        } else {
            output <- input %>% filter(comp != "Ref") %>% 
                gather(size_var,value,area,ampl) %>% 
                group_by(row,comp,class) %>% 
                mutate(acceptable_peak_value = value[which(size_var == acceptable_peak_units)]) %>% 
                filter(acceptable_peak_value > smallest_acceptable_peak & acceptable_peak_value < largest_acceptable_peak) %>% 
                mutate(size_value = value[which(size_var == size_filter)],
                       size_lower = NA,
                       size_upper = NA,
                       size_group = NA) %>% 
                spread(size_var,value) %>% 
                mutate(dD_presize = dD_processing,
                       dD_processing = dD_presize,
                       dD_size = dD_processing) %>% 
                filter(!is.na(dD_processing))
            
            list("output" = output,
                 "size_filter" = NULL,
                 "size_effect_raw" = NULL,
                 "size_raw_bycomp_plot" = ggplot()+annotate("text",x=0,y=0,label="No Size Effect Standards Found.",color="red")+theme_void(),
                 "size_raw_grouped_plot" = ggplot()+annotate("text",x=0,y=0,label="No Size Effect Standards Found.",color="red")+theme_void(),
                 "size_corrected_bycomp_plot" = ggplot()+annotate("text",x=0,y=0,label="No Size Effect Standards Found.",color="red")+theme_void(),
                 "size_corrected_grouped_plot" = ggplot()+annotate("text",x=0,y=0,label="No Size Effect Standards Found.",color="red")+theme_void(),
                 "size_correction_plot" = ggplot()+annotate("text",x=0,y=0,label="No Size Effect Standards Found.",color="red")+theme_void()
            )
        }
}

normalization_function <- function(input,normalization_option,normalization_comps) {
    
    scale_check <- input %>% filter(grepl("standard",id2))
    
    if(length(scale_check$id2) > 0){
        normalization_option = ifelse(grepl("Linear interpolation between adjacent normalization standards",normalization_option),1,2)
        mix_to_use = input %>% filter(grepl("standard",id2)) %>% ungroup() %>% select(std_mix) %>% distinct() %>% .$std_mix
        comps_to_use = input %>% filter(std_mix == mix_to_use & comp_class %in% normalization_comps) %>% ungroup() %>% select(comp_class) %>% distinct() %>% .$comp_class
        
        scale_normalization_1 <- input %>% filter(comp != "Ref") %>% 
            filter(grepl("standard",id2) & comp_class %in% comps_to_use) %>% 
            filter(!is.na(dD_known)) %>% 
            select(row,row_base,std_mix,id2,comp,class,dD_known,dD_processing) %>% 
            group_by(comp,class) %>% 
            mutate(i = dense_rank((row_base))) %>%  # Generate an indexing variable for each set of standard injections... only relevant using normalization_option 1
            group_by(i) %>%
            # The next mutate is doing the heavy lifting. Each 'bin' is a given drift sample (i) and the next drift sample (i+1)
            # We are basically just doing a two point linear regression to find the slope and y-intercept of the line.
            mutate(slope_1 = lm(dD_known~dD_processing,data=.[which(.$i == i[1] | .$i == i[1]+1 ),])$coefficients[2],
                   intercept_1 = lm(dD_known~dD_processing,data=.[which(.$i == i[1] | .$i == i[1]+1 ),])$coefficients[1]) %>% 
            ungroup() %>% 
            mutate(slope_1 = ifelse(i == max(i),slope_1[which(i == max(i)-1)],slope_1), # Assigns second to last slope to the last drift set.
                   intercept_1 = ifelse(i == max(i),intercept_1[which(i == max(i)-1)],intercept_1), # Same as above.
                   slope_2 = lm(dD_known~dD_processing)$coefficients[2], # Calculates the slope for normalization_option = 2
                   intercept_2 = lm(dD_known~dD_processing)$coefficients[1])  # Calculates the intercept for normalization_option = 2
        
        # Next, we need to make a data frame to fill in all the rows (injections) that occurred between the drift samples.
        scale_normalization_2 <- scale_normalization_1 %>% 
            full_join(data.frame(row = seq(1,max(input$row))),by="row") %>% 
            gather(coef,value,slope_1,slope_2,intercept_1,intercept_2) %>% 
            separate(coef,c("coef","normalization_method"),sep="_") %>% 
            select(row,coef,normalization_method,value) %>% distinct() %>%
            filter(normalization_method %in% normalization_option) %>%
            group_by(coef) %>% 
            mutate(value = ifelse(row == min(row),value[which(row == min(row[which(!is.na(value))]))],value)) %>% 
            arrange(coef,row) %>%  # Sort the dataset by coefficient and then by row. This is important because we use a 'last-observation-carried-forward' command...
            mutate(value = na.locf(value)) %>% # ... which fills each row with the preceding row's coefficient value.
            distinct() %>% 
            spread(coef,value)
        
        normalization_plot <- scale_normalization_1 %>% 
            gather(coef,value,slope_1:intercept_2) %>% 
            separate(coef,c("coef","normalization_method"),sep="_") %>% 
            spread(coef,value) %>% 
            group_by(i,normalization_method) %>% 
            mutate(method_label = factor(normalization_method,levels=c(1,2),labels=c("Bracketed Normalization","Full-Run Normalization")),
                   equation = paste0("y = ",round(slope[1],4)," * x + ",round(intercept[1],4)),
                   label = paste(method_label,equation,sep="\n")) %>% 
            filter(normalization_method %in% normalization_option) %>% 
            ggplot(aes(x=dD_processing, y=dD_known)) +
            geom_point() +
            geom_smooth(method="lm",se=F) +
            facet_rep_wrap(~label) +
            labs(x="Observed \u04B4D ( \u2030 )",
                 y="Known \u04B4D ( \u2030 )",
                 title="VSMOW-SLAP Scale Normalization Curve")
        
        output <- input %>% filter(comp != "Ref") %>% 
            full_join(scale_normalization_2,by = "row") %>% 
            mutate(dD_prenormalization = dD_processing,
                   dD_processing = dD_prenormalization * slope + intercept,
                   dD_normalization = dD_processing) %>% 
            anti_join(scale_normalization_1,by = c("row", "id2", "comp", "class", "row_base", "std_mix", "dD_known"))
        
        list("output" = output,
             "scale_normalization_1" = scale_normalization_1,
             "scale_normalization_2" = scale_normalization_2,
             "normalization_plot" = normalization_plot)
    } else {
        output <- input %>% filter(comp != "Ref") %>% 
            mutate(dD_prenormalization = dD_processing,
                   dD_processing = dD_prenormalization,
                   dD_normalization = dD_processing)
        list("output" = output,
             "scale_normalization_1" = NULL,
             "scale_normalization_2" = NULL,
             "normalization_plot" = ggplot()+annotate("text",x=0,y=0,label="No Normalization Standards Found!",color="red")+theme_void())
    }
}

control_function <- function(input,normalization_comps,processing_order) {
    control_check <- input %>% filter(grepl("control",id2))
    
    if(length(control_check$id2) > 0) {
        control_standards <- input %>% 
            filter(grepl("control",id2)) %>% 
            filter(!(comp_class %in% normalization_comps)) %>% 
            ungroup() %>% 
            select(comp_class,conc,std_mix,dD_processing,dD_known) %>% 
            mutate(ungrouped_full_rmse = round(sqrt(sum((dD_known - dD_processing)^2)/(n()-1)),1)) %>% 
            group_by(std_mix,conc) %>% 
            mutate(std_mix_rmse = round(sqrt(sum((dD_known - dD_processing)^2)/(n()-1)),1)) %>% 
            group_by(comp_class,conc,std_mix,ungrouped_full_rmse,std_mix_rmse) %>%
            summarize(observed_value = round(mean(dD_processing),1),
                      observed_sd = round(sd(dD_processing),1),
                      accepted_value = round(unique(dD_known)[1],1),
                      root_mean_square_error = round(sqrt(sum((dD_known - dD_processing)^2)/(n()-1)),1),
                      mean_signed_difference = round(sum(dD_known-dD_processing)/(n()-1),1)) %>% 
            arrange(std_mix,conc,comp_class)
        
        proc_table <- data.frame(levels = c("dD_raw","dD_drift","dD_size","dD_normalization"),
                                 proc_label = c("Raw","Drift","Size","Scale Normalization")) %>% 
            full_join(data.frame(proc_label = c("Raw",processing_order),
                                 order = 1:4),by = "proc_label") %>% 
            mutate(label = paste0("Step ",order,": ",proc_label)) %>% 
            arrange(order)
        
        correction_visual_plot <- input %>% 
            filter(!is.na(dD_known)) %>% 
            filter(!grepl("derivatization",id2)) %>% 
            gather(dD_type,value,dD_raw,dD_drift,dD_size,dD_normalization) %>% 
            mutate(dD_type = factor(dD_type,levels=proc_table$levels,ordered=T,
                                    labels=proc_table$label)) %>% 
            group_by(dD_type,dD_known,dD_known_sd,std_mix,conc,comp) %>% 
            summarize(dD_mean = mean(value),
                      dD_sd = sd(value)) %>% 
            mutate(std_conc = paste(std_mix,conc))
        
        corrections_plot <- ggplot(correction_visual_plot, aes(x=dD_mean,y=dD_known)) +
            geom_errorbar(aes(ymin = dD_known - dD_known_sd, ymax = dD_known + dD_known_sd),width=0) +
            geom_errorbarh(aes(xmin = dD_mean - dD_sd,xmax = dD_mean + dD_sd),height=0) +
            geom_smooth(method="lm",se=F,color="black",size=0.25) +
            geom_point(aes(color=std_conc),size=2.5) +
            facet_rep_wrap(~dD_type) + 
            labs(x = "\u04B4D Observed ( \u2030 )",
                 y = "\u04B4D Known ( \u2030 )",
                 color = "Standard + Concentration")
        
        list("control_standards" = control_standards,
             "corrections_plot" = corrections_plot)
    } else {
        list("control_standards" = NULL,
             "corrections_plot" = ggplot()+annotate("text",x=0,y=0,label="No Control (i.e., QA/QC) Standards Found.",color="red")+theme_void())
    }
        
}

derivatization_correction <- function(input,derivatization_table,derivatization_option) {
    
    if(!is.null(input)){
        derivative_check <- input %>% filter(grepl("derivatization",id2))
        
        derivatization_choice = ifelse(derivatization_option == "Template-defined derivative \u03B4D.","template",
                                ifelse(derivatization_option == "Derivatization standard in sequence \u03B4D.","inrun",
                                ifelse(derivatization_option == "Do not correct for derivative hydrogen.","disabled",
                                       "derivatization_option in derivatization_correction broken! Check that the radioButton choices still match!")))
        
        if(length(derivative_check$id2) > 0) {
            inrun <- derivative_check %>% left_join(derivatization_table[,1:4],by = c("comp","class")) %>% 
                ungroup() %>% 
                mutate(total_hydrogen_count = bound_hydrogen_count + derivative_hydrogen_count,
                       derivative_dD_inrun = (total_hydrogen_count*dD_processing - bound_hydrogen_count*dD_known) / derivative_hydrogen_count,
                       error_counts = ifelse(anyNA(total_hydrogen_count),"[ No match found for derivatization standard in Derivatization Table sheet! ]",""),
                       error_known = ifelse(anyNA(dD_known),"[ No match found for derivatization standard in Standards sheet! ]",""),
                       error_known_sd = ifelse(anyNA(dD_known_sd) & is.null(error_known),"[ Missing uncertainty of known value for derivatization standard in Standards sheet!",""),
                       errors = paste0(error_counts,error_known,error_known_sd),
                       derivatization_filter = "inrun") %>% 
                group_by(derivatization_filter,errors) %>% 
                summarize(derivative_dD_inrun_mean = round(mean(derivative_dD_inrun),1),
                          derivative_dD_inrun_precision = sd(derivative_dD_inrun),
                          derivative_dD_inrun_knownerror = mean(dD_known_sd),
                          derivative_dD_inrun_uncertainty = round(sqrt(derivative_dD_inrun_precision^2 + derivative_dD_inrun_knownerror^2),1)) %>% 
                rename(derivative_dD = derivative_dD_inrun_mean,
                       derivative_dD_uncertainty = derivative_dD_inrun_uncertainty) %>% 
                select(derivative_dD,derivative_dD_uncertainty,derivatization_filter,errors)
        } else {
            inrun <- data.frame(derivative_dD = NA,
                                derivative_dD_uncertainty = NA,
                                derivatization_filter = "inrun",
                                errors = "[ No derivatization standards identified in Identifier 2! ]")
        }
        
        derivatives_to_use <- derivatization_table[1,c(5:6)] %>% 
            mutate(derivatization_filter = "template",
                   derivative_dD_error = ifelse(is.na(derivative_dD),"[ Missing derivative_dD in Derivatization sheet! ]",""),
                   derivative_dD_uncertainty_error = ifelse(is.na(derivative_dD_uncertainty),"[ Missing derivative_dD_uncertainty in Derivatization sheet! ]",""),
                   errors = paste0(derivative_dD_error,derivative_dD_uncertainty_error)) %>% 
            select(-c(derivative_dD_error,derivative_dD_uncertainty_error)) %>% 
            rbind(inrun) %>% 
            rbind(data.frame(derivative_dD = NA, derivative_dD_uncertainty = NA, derivatization_filter = "disabled",errors = "")) %>% 
            filter(derivatization_filter %in% derivatization_choice) %>% 
            mutate(derivatization_filter = derivatization_option)
        
        sample_results <- left_join(input,derivatization_table[,1:4],by = c("comp", "class")) %>% 
            group_by(row,comp,class) %>% 
            mutate(derivative_dD = as.numeric(derivatives_to_use$derivative_dD),
                   derivative_dD_uncertainty = as.numeric(derivatives_to_use$derivative_dD_uncertainty)) %>% 
            mutate(bound_hydrogen_count = ifelse(is.na(bound_hydrogen_count),1,bound_hydrogen_count),
                   derivative_hydrogen_count = ifelse(is.na(derivative_hydrogen_count),0,derivative_hydrogen_count),
                   total_hydrogen_count = bound_hydrogen_count + derivative_hydrogen_count,
                   derivative_dD_uncertainty = ifelse(derivative_hydrogen_count == 0,0,derivative_dD_uncertainty),
                   final_dD = (total_hydrogen_count*dD_processing - derivative_hydrogen_count*derivative_dD)/bound_hydrogen_count,
                   final_dD = ifelse(derivatization_choice == "disabled",dD_processing,final_dD),
                   final_dD = ifelse(id2 == "sample",final_dD,dD_processing))
        
        list("output" = sample_results,
             "derivative_dD_selected" = derivatives_to_use)
    } else {
        list("output" = NULL,
             "derivative_dD_selected" = NULL)
    }
}

final_sample_function <- function(input,control_error,control_option,control_standards) {
    
    
}

# UI
{
    ui <- dashboardPage(
        skin = "green",
        dashboardHeader(title = "GC-IRMS \u03B4D"),
        dashboardSidebar(
            sidebarMenu(
                menuItem("Ingest Data",tabName = "ingest_tab",icon = icon("table")),
                menuItem("Peaks & Diagnostics",tabName = "peak_tab",icon = icon("chart-line")),
                menuItem("Processing Order",tabName = "order_tab",icon = icon("list-ol")),
                menuItem("Drift Correction",tabName = "drift_tab",icon = icon("shoe-prints")),
                menuItem("Size Correction",tabName = "size_tab",icon = icon("signal")),
                menuItem("Scale Normalization",tabName = "normalization_tab", icon = icon("balance-scale")),
                menuItem("Control Standards",tabName = "control_tab", icon = icon("bullseye")),
                menuItem("Derivatization",tabName = "derivatization_tab", icon = icon("vials")),
                menuItem("Sample Error",tabName = "error_tab", icon = icon("square-root-alt"))
            )
        ),
        dashboardBody(
            tags$style("* {font-family: Arial;}"),
            tags$h1(tags$style("h1 {font-family: Arial;}")),
            tags$h3(tags$style("h3 {font-family: Arial;}")),
            tabItems(
                tabItem(tabName = "ingest_tab",
                        fluidPage(
                            box(title = h1("Shiny Post-Processor: Compound Specific \u03B4D via GC-HTC-IRMS", strong("(expand for instructions!)"),hr()),
                                solidHeader=T,
                                width=12,
                                collapsible = T,
                                collapsed = T,
                                uiOutput("introduction")),
                            box(title = "Select IRMS export file:",
                                width=3,
                                fileInput("raw_irms_file",
                                          label = "May be a CSV, XLS, or XLSX, but headers must be identified in the template if not an Isodat 3.0 CSV export.",
                                          accept = c(".csv",".xlsx",".xls"))),
                            box(title = "Select GC-IRMS Template file:",
                                width=3,
                                fileInput("gcirms_template",
                                          label = "This must follow the structure of the original template! Do not delete any sheets!",
                                          accept = c(".xlsx"))),
                            column(3,
                                   box(width=12,height=85,actionButton("load_demos",label="Loads examples from Github.",icon = icon("file"))),
                                   infoBoxOutput("gcirms_template_status",width=12)),
                            column(3,
                                infoBoxOutput("raw_irms_status",width=12),
                                infoBoxOutput("raw_irms_check",width=12)),
                            box(title = "Raw IRMS Export",
                                width=12,
                                collapsible = T,
                                collapsed = T,
                                DTOutput("raw_irms_data"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll"),
                            box(title = "Sample Info",
                                width=12,
                                collapsible = T,
                                collapsed = T,
                                DTOutput("sample_info"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll"),
                            box(title = "Standard Info",
                                width=12,
                                collapsible = T,
                                collapsed = T,
                                DTOutput("standard_info"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll"),
                            box(title = "Initial IRMS Clean-up",
                                width=12,
                                collapsible = T,
                                collapsed = T,
                                DTOutput("ingestdata"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll")
                        )
                ),
                tabItem(tabName = "peak_tab",
                        fluidPage(
                            box(title = NULL,
                                width=6,
                                selectInput("compound_option",
                                            label = "Select how to assign compound names to peaks:",
                                            choices = c("Assigned by IRMS export in Component/comp column.",
                                                        "Assigned by retention time matching with Retention Times sheet."))),
                            box(title="Retention Times",
                                width=6,
                                collapsible = T,
                                collapsed = T,
                                DTOutput("rt_table"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll"),
                            box(title="Uncorrected IRMS Data with Peak Assignments",
                                width=12,
                                collapsible = T,
                                collapsed = T,
                                DTOutput("data"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll"),
                            box(title = "IRMS Drift Based on Reference Peaks",
                                width=4,
                                plotOutput("plot_irms_drift")),
                            box(title = "Peak Area to Injection Concentration",
                                width=4,
                                plotOutput("plot_area_curve")),
                            box(title = "Peak Area to Injection Concentration",
                                width=4,
                                plotOutput("plot_ampl_curve"))
                        )
                ),
                tabItem(tabName = "order_tab",
                        fluidPage(
                            box(title = "Select the order to perform corrections. Each correction can also be disabled by an option in its tab.",
                                width=12),
                            box(title = NULL,
                                width=4,
                                selectInput("first_correction",
                                            label = "First Correction",
                                            choices = c("Drift","Size","Scale Normalization"))),
                            box(title = NULL,
                                width=4,
                                selectInput("second_correction",
                                            label = "Second Correction",
                                            choices = c("Drift","Size","Scale Normalization"))),
                            box(title = NULL,
                                width=4,
                                selectInput("third_correction",
                                            label = "Third Correction",
                                            choices = c("Drift","Size","Scale Normalization")))
                        )
                ),
                tabItem(tabName = "drift_tab",
                        fluidPage(
                            fluidRow(
                                box(title = NULL,
                                    width = 6,
                                    radioButtons("drift_option",
                                                 label = "Select the drift type to apply:",
                                                 choices = c("Linear interpolation between adjacent drift samples (use this when drift appears non-linear)",
                                                             "Linear regression across all drift samples (use this when drift appears linear)",
                                                             "No drift correction (use this when there is no apparent drift or if you use bracketed scale normalization)"))),
                                box(title = NULL,
                                    width = 6,
                                    radioButtons("drift_comp",
                                              label = "Select the compound to use for drift correction:",
                                              choices = "Mean drift of all compounds"))
                            ),
                            box(title = NULL,
                                width = 4,
                                plotOutput("drift_raw_plot")),
                            box(title = NULL,
                                width = 4,
                                plotOutput("drift_corrected_standards_plot")),
                            box(title = NULL,
                                width = 4,
                                plotOutput("drift_correction_by_row_plot")),
                            box(title = "Drift Correction Output",
                                width=12,
                                DTOutput("drift_corrected_data"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll")
                        )
                ),
                tabItem(tabName = "size_tab",
                        fluidPage(
                            fluidRow(
                                box(title = NULL,
                                    width = 6,
                                    radioButtons("size_option",
                                                 label = "Select the size effect independent variable to use:",
                                                 choices = c("Peak Height (amplitude, mV)",
                                                             "Peak Area (Vs)",
                                                             "No size effect correction.")),
                                    numericInput("size_cutoff",
                                              label = "Enter a cut-off value to perform a segmented size effect correction. Leave blank to use a single regression",
                                              value = NA,
                                              step = 10),
                                    radioButtons("size_small_peak_action",
                                                 label = "Select what action to take for 'Small' size peaks.",
                                                 choices = c("Remove observed size effect.",
                                                             "No size effect correction.")),
                                    radioButtons("size_normal_peak_action",
                                                 label = "Select what action to take for 'Normal' size peaks.",
                                                 choices = c("Remove observed size effect.",
                                                             "No size effect correction."))),
                                box(title = NULL,
                                    width = 6,
                                    radioButtons("acceptable_peak_units",label = "Select the unit for the smallest/largest peak filters.",
                                                 choices = c("Peak Height (amplitude, mV)",
                                                             "Peak Area (Vs)")),
                                    numericInput("smallest_acceptable_peak",
                                                 label = "Enter a value to remove any peak smaller than this. Leave blank to disable lower peak limit.",
                                                 value = NA),
                                    numericInput("largest_acceptable_peak",
                                                 label = "Enter a value to remove any peak larger than this. Leave blank to disable upper peak limit.",
                                                 value = NA)),
                                box(title = NULL,
                                    width = 3,
                                    radioButtons("size_toosmall_peak_action",
                                                 label = "Select what action to take for peaks smaller than the lower range of the size effect curve.",
                                                 choices = c("Use 'Small' size effect function.",
                                                             "No size effect correction.",
                                                             "Remove these from results."))),
                                box(title = NULL,
                                    width = 3,
                                    radioButtons("size_large_peak_action",
                                                 label = "Select what action to take for peaks larger than the upper range of the size effect curve.",
                                                 choices = c("Remove size effect with 'Normal' size effect function.",
                                                             "No size effect correction.",
                                                             "Remove these from results.")))
                            ),
                            box(title = NULL,
                                width = 6,
                                plotOutput("size_raw_grouped_plot")),
                            box(title = NULL,
                                width = 6,
                                plotOutput("size_corrected_grouped_plot")),
                            box(title = NULL,
                                width = 12,
                                plotOutput("size_correction_plot")),
                            box(title = NULL,
                                width = 12,
                                plotOutput("size_raw_bycomp_plot")),
                            box(title = NULL,
                                width = 12,
                                plotOutput("size_corrected_bycomp_plot")),
                            box(title = "Size Correction Output",
                                width=12,
                                DTOutput("size_corrected_data"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll")
                        )
                 ),
                 tabItem(tabName = "normalization_tab",
                        fluidPage(
                            fluidRow(
                                box(title = NULL,
                                    width = 6,
                                    radioButtons("normalization_option",
                                                 label = "Select the desired method to perform scale normalization:",
                                                 choices = c("Linear interpolation between adjacent normalization standards (use this if drift-correction is untenable)",
                                                             "Linear regression across all normalization standards (use this if you drift-correct)"))),
                                box(title = NULL,
                                    width = 6,
                                    checkboxGroupInput("normalization_comps",
                                                       label = "Check at least 2 compounds to use as basis for scale normalization. If not defined in the Initials tab, the two compounds that bracket all the compounds in the scale normalization standard are used."))
                            ),
                            fluidRow(title = NULL,
                                column(width=6,offset=3,align="center",plotOutput("normalization_plot",height = 600))),
                            box(title = "Scale Normalized Output",
                                width=12,
                                DTOutput("scale_normalization_corrected_data"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll")
                        )
                ),
                tabItem(tabName = "control_tab",
                        fluidPage(
                            box(title = NULL,
                                width = 12,
                                DTOutput("control_standards")),
                            fluidRow(title = NULL,
                                     column(width=6,offset=3,align="center",plotOutput("corrections_plot",height = 600)))
                        )
                ),
                tabItem(tabName = "derivatization_tab",
                        fluidPage(
                            box(title = NULL,
                                width = 3,
                                radioButtons("derivatization_option",
                                             label = "Select the desired source of the derivative hydrogen \u03B4D:",
                                             choices = c("Template-defined derivative \u03B4D.",
                                                         "Derivatization standard in sequence \u03B4D.",
                                                         "Do not correct for derivative hydrogen."))),
                            box(title = "Derivative Hydrogen Table",
                                width = 9,
                                collapsible = T,
                                collapsed = T,
                                DTOutput("derivatization_table"),
                                style="height:250px; overflow-y: scroll;overflow-x: scroll"),
                            box(title = "Selected Derivative \u03B4D Values",
                                width = 12,
                                tableOutput("derivative_dD_selected")),
                            box(title = "Sample Derivatization Correction",
                                width = 12,
                                DTOutput("sample_derivatization_table"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll")
                        )
                ),
                tabItem(tabName = "error_tab",
                        fluidPage(
                            box(title = NULL,
                                width = 6,
                                radioButtons("accuracy_source",
                                             label = "Choose which pool of control standards to use for accuracy propagation:",
                                             choices = c("Sequence control standards only.",
                                                         "Long-term control standards only.",
                                                         "Sequence & long-term control standards."))),
                            box(title = NULL,
                                width = 6,
                                checkboxGroupInput("accuracy_standards",
                                                   label = "Choose which standards to calculate RMS error from:")),
                            box(title = "Chosen Error Values & Warnings",
                                width = 12,
                                tableOutput("chosen_error_table")),
                            box(title = "Final Sample Errors",
                                width = 12,
                                DTOutput("final_sample_error_table"),
                                style="height:500px; overflow-y: scroll;overflow-x: scroll")
                        )
                )
            )
        )
    )
}

server <- function(input, output, session) {
   #Ingest Data Tab
    {
        output$introduction <- renderUI({
            tagList(
                p("Source code, template, and example data can be found ",a("on Github.",href="https://github.com/JackAHutchings/Shiny-Post-Processor")," For development questions/bug
                  reports contact the author, Jack Hutchings, via Github."),
                h3("Sequence Construction"),
                p("When setting up your sample sequence, you *must* reserve two columns for 'Identifier 1' and 'Identifier 2'. A third, optional column, is 'Preparation'.",
                  style = "padding-left:2em"),
                p(strong("Identifier 1:")," This is the primary identifier of the vial you are injecting. If the vial is a standard, then there should be matching entries in your
                  template file (Standards sheet; see below) to fill in known values. If the vial is a sample, then this can be whatever value you desire. However, if you are
                  using the optional Samples sheet in the template file, then you'll want these to match so you have proper metadata.",
                  style = "padding-left:2em"),
                p(strong("Identifier 2:")," This is reserved (i.e., should only contain values as stated here) and must be used to describe the role of the injection. The 
                  following roles are understood:",
                  style = "padding-left:2em"),
                tags$ul(
                    tags$li("sample (unknown)",style = "list-style-position: inside;text-indent: -1em;padding-left: 1em;"),
                    tags$li("standard (scale normalization standard)",style = "list-style-position: inside;text-indent: -1em;padding-left: 1em;"),
                    tags$li("control (known standard to estimate accuracy)",style = "list-style-position: inside;text-indent: -1em;padding-left: 1em;"),
                    tags$li("drift (account for instrument drift)",style = "list-style-position: inside;text-indent: -1em;padding-left: 1em;"),
                    tags$li("size (peak size to isotope ratio correction)",style = "list-style-position: inside;text-indent: -1em;padding-left: 1em;"),
                    tags$li("derivatization (calculate added hydrogen during derivatization)",style = "list-style-position: inside;text-indent: -1em;padding-left: 1em;"),
                    tags$li("warm-up (discarded injections)",style = "list-style-position: inside;text-indent: -1em;padding-left: 1em;")
                    ),
                p("Standards may serve multiple roles (i.e, drift and control) that should be separated by an underscore (i.e., drift_control) in Identifier 2.
                  If a standard is assigned both the standard and control role (i.e., standard_control), then the compounds selected for scale normalization will be
                  excluded from use as control standards. Only sample and standard roles are required to be used. If you have pre-determined your derivative hydrogen, 
                  the known value can be entered in the Excel template file (see below)",
                  style = "padding-left:2em"),
                p(strong("Preparation:")," (optional) This column is used to indicate the concentration ( ng/\u03BCL ) of a vial. This is only necessary if size correction is desired 
                  or if multiple concentrations of one of the other standards are used. Sample concentration is not interpreted.",
                  style = "padding-left:2em"),
                h3("GCIRMS Template Completion"),
                p("The Microsoft Excel template file can be found at",a("this app's Github repository.",href="https://github.com/JackAHutchings/Shiny-Post-Processor/tree/master/GCIRMS_Processor"),
                  " The file contains instructions (spreadsheet cells with bright yelow backgrounds), but general usage of each 'Sheet' within the file is described below. Each sheet is labeled
                  as (Optional) or (Required) if it is necessary for the user to change values in the sheet for functioning. However, even if a sheet is optional,",
                  strong("do not delete or rename any sheets!")," You may add sheets as desired, but all the sheets and their names in the original template must be present for the app 
                  to function.",
                  style = "padding-left:2em"),
                tags$ul(
                    tags$li(strong("Sequence Table:")," (Optional) The sequence table is designed specifically for copy/pasting into an Isodat 3.0 GC Isolink II Sequence file. If you are running 
                      an older version or some other software, this template will not be particularly useful. However, the app does not use this sheet, so you may delete it or replace it with 
                      a template of your own.",
                            style = "list-style-position: inside;text-indent: -20px;padding-left: 20px;"),
                    tags$li(strong("Samples:")," (Optional) If used, the 'Identifier 1' column should match 'Identifier 1' in your sequence table and IRMS output. The remaining
                      columns (B through G) can have their headers changed and contents filled with whatever sample metadata you require.",
                            style = "list-style-position: inside;text-indent: -20px;padding-left: 20px;"),
                    tags$li(strong("Standards:"), " (Required) Standards with known compositions (scale normalization and controls) must have their known values entered here. You may enter
                      as many standards as you desire; only those with matching Identifier 1 values will be used by the script. If a standard contains compounds at different concentrations,
                      then their concentation relative to 'Preparation' in your Sequence Table/IRMS Export should be entered here. If you don't have or want to use a 'Preparation' column, the actual
                      concentrations could be entered here. To have a standard at multiple concentrations (e.g., your peak size effect standard), you would need each concentration level to have its
                      own Identifier 1 and concentration information entered in relative_concentration.",
                            style = "list-style-position: inside;text-indent: -20px;padding-left: 20px;"),
                    tags$li(strong("Headers:")," (Required) This sheet informs the app what the actual header values are in your IRMS export file. Some of these are optional, but for full functionality it is
                      recommended to include as many as you have available. If you are exporting using Isodat 3.0 as a CSV, then you only need to select the correct headers in your export template,
                      as the values here match those. If you are using older/other software, you will likely need to identify each column and provide the correct header in your IRMS export file.",
                            style = "list-style-position: inside;text-indent: -20px;padding-left: 20px;"),
                    tags$li(strong("Retention Times:")," (Optional) A list of compound names & classes and their GC retention times. If you use your software's built-in RT table and your compounds are
                      already identified in your IRMS export, then this is not needed.",
                            style = "list-style-position: inside;text-indent: -20px;padding-left: 20px;"),
                    tags$li(strong("Derivatization:")," (Optional) A list of compound names & classes and their hydrogen atom counts. If you are running non-derivatized compounds (e.g., n-alkanes), then
                      this table is not necessary. The known value of your derivative \u03B4D should be entered here with appropriate error. This may be left empty if you include a derivatization
                            standard with the sequence.",
                            style = "list-style-position: inside;text-indent: -20px;padding-left: 20px;"),
                    tags$li(strong("Initials:")," (Optional) The initial choices for each processing option. If you have a preferred method, then you may want to
                      change these to match your preference. To work, these must be exactly the same as the options visible throughout the app.",
                            style = "list-style-position: inside;text-indent: -40px;padding-left: 40px;")
                    ),
                h3("IRMS Export Formatting"),
                p("Your IRMS export should be either a CSV, XLS, or XLSX file. The required columns in the export are identified in the Headers sheet of the Excel template. You must have
                  these headers present ",strong("and")," they must match the values in the 'export_header' column of the Headers sheet.",
                  style = "padding-left:2em"),
                h3("Usage & Export"),
                p("Once you have successfully imported your IRMS data and an appropriately completed template, simply progress through the app using the sidebar and make the desired choices.
                  Export options are still in progress!",
                  style = "padding-left:2em")
            )
        })
        
        # demo_paths <- eventReactive(input$load_demos,{
        #            c("https://raw.githubusercontent.com/JackAHutchings/Shiny-Post-Processor/master/GCIRMS_Processor/Example%20Export.csv",
        #              "https://github.com/JackAHutchings/Shiny-Post-Processor/blob/master/GCIRMS_Processor/GCIRMS%20Template.xlsx?raw=true")
        # })
        
        
        ingest <- reactive({
            
            # if ( !is.null(demo_paths()) ) {
            #     GET(demo_paths()[2],write_disk(tf <- tempfile(fileext = ".xlsx")))
            #     ingest_function(demo_paths()[1],tf)}
            if ( is.null(input$raw_irms_file) | is.null(input$gcirms_template) ) {return(NULL)} else{
            if (!(is.null(input$raw_irms_file) & is.null(input$gcirms_template)) &
                    (grepl("csv",input$raw_irms_file$datapath) | grepl("xls",input$raw_irms_file$datapath)) &
                    grepl("xlsx",input$gcirms_template$datapath)){ingest_function(input$raw_irms_file$datapath,input$gcirms_template$datapath)}
            }
            
        })
        
        cal_comp_list <- reactive({
            if(!is.null(ingest()$output)){data()$data %>% filter(comp != "Ref") %>%  # Removes the reference peaks.
                    filter(grepl("standard",id2)) %>%
                    select(comp_class,dD_known) %>% distinct()} else return(NULL)})
        default_cal_comp <- reactive({
            if(!is.null(cal_comp_list())){
                if(anyNA(ingest()$initials$normalization_comps)){cal_comp_list() %>% filter(dD_known == max(dD_known) | dD_known == min(dD_known)) %>% .$comp_class}
                if(!anyNA(ingest()$initials$normalization_comps)){ingest()$initials$normalization_comps}}else return(NULL)})

        observe({
            updateSelectInput(session,"first_correction",selected = ingest()$initials$first_correction[1])
            updateSelectInput(session,"second_correction",selected = ingest()$initials$second_correction[1])
            updateSelectInput(session,"third_correction",selected = ingest()$initials$third_correction[1])
            updateRadioButtons(session,"drift_option",selected = ingest()$initials$drift_option[1])
            updateRadioButtons(session,"drift_comp",selected = ingest()$initials$drift_comp[1])
            updateRadioButtons(session,"size_option",selected = ingest()$initials$size_option[1])
            updateNumericInput(session,"size_cutoff",value = ingest()$initials$size_cutoff[1])
            updateRadioButtons(session,"size_toosmall_peak_action",selected = ingest()$initials$size_toosmall_peak_action[1])
            updateRadioButtons(session,"size_small_peak_action",selected = ingest()$initials$size_small_peak_action[1])
            updateRadioButtons(session,"size_normal_peak_action",selected = ingest()$initials$size_normal_peak_action[1])
            updateRadioButtons(session,"size_large_peak_action",selected = ingest()$initials$size_large_peak_action[1])
            updateRadioButtons(session,"acceptable_peak_units",selected = ingest()$initials$acceptable_peak_units[1])
            updateNumericInput(session,"largest_acceptable_peak",value = ingest()$initials$largest_acceptable_peak[1])
            updateNumericInput(session,"smallest_acceptable_peak",value = ingest()$initials$smallest_acceptable_peak[1])
            updateRadioButtons(session,"normalization_option",selected = ingest()$initials$normalization_option[1])
            updateCheckboxGroupInput(session,"normalization_comps",choices = cal_comp_list()$comp_class, selected = default_cal_comp())
            updateRadioButtons(session,"derivatization_option",selected = ingest()$initials$derivatization_option[1])
            },priority = 1)
        observe({updateSelectInput(session,"compound_option",selected = ingest()$initials$compound_option[1])})
        
        output$raw_irms_status <- renderInfoBox({
            if (is.null(input$raw_irms_file)) {text = "No File Uploaded!"; use_color = "blue"}
            else if( !(grepl("csv",input$raw_irms_file$datapath)|grepl("xls",input$raw_irms_file$datapath) )) {text = "Raw IRMS file must be CSV, XLS, or XLSX!"; use_color = "red"}
            else if( !is.null(input$raw_irms_file) & (grepl("csv",input$raw_irms_file$datapath)|grepl("xls",input$raw_irms_file$datapath))) {text = "File Uploaded."; use_color = "green"}
            
            infoBox("Status:",text,icon = icon("file-alt"),color = use_color)
        })
        output$gcirms_template_status <- renderInfoBox({
            if (is.null(input$gcirms_template)) {text = "No File Uploaded!"; use_color = "blue"}
            else if( !grepl("xlsx",input$gcirms_template$datapath) ) {text = "GCIRMS Template file must be XLSX!"; use_color = "red"}
            else if (!is.null(input$gcirms_template) & grepl("xlsx",input$gcirms_template$datapath)) {text = "File Uploaded."; use_color = "green"}
            infoBox("Status:",text,icon = icon("file-excel"),color = use_color)
        })
        
        output$raw_irms_check <- renderInfoBox({
            if(is.null(ingest())){use_text = "Waiting for both file uploads..."; use_color = "blue"; use_icon = icon("ellipsis-h")}
            if(isTRUE(ingest()$raw_irms_check)){use_text = "Inital File Check Passed."; use_color = "green"; use_icon = icon("check")}
            if(is.character(ingest()$raw_irms_check)){use_text = paste("Header Mismatch: ",paste(ingest()$raw_irms_check,collapse=", ")); use_color = "red"; use_icon = icon("exclamation-triangle")}
            infoBox("File Check:",value=use_text,icon=use_icon,color=use_color)
        })
        
        output$raw_irms_data <- renderDT(ingest()$raw_irms_data,options=list('lengthMenu'=JS('[[10,25,50,-1],[10,25,50,"All"]]'),searching=FALSE),class='white-space:nowrap')
        output$sample_info <- renderDT(ingest()$sample_info,options=list('lengthMenu'=JS('[[10,25,50,-1],[10,25,50,"All"]]'),searching=FALSE),class='white-space:nowrap')
        output$standard_info <- renderDT(ingest()$standard_info,options=list('lengthMenu'=JS('[[10,25,50,-1],[10,25,50,"All"]]'),searching=FALSE),class='white-space: nowrap')
        output$ingestdata <- renderDT({ingest()$output},options=list('lengthMenu'=JS('[[10,25,50,-1],[10,25,50,"All"]]'),searching=FALSE),class='white-space:nowrap')
        
    }
   # Compound Assignment & Initial Diagnostic Plots Tab
    {
        data <- reactive({if(is.null(ingest()$output)){NULL}else{comp_assign_function(ingest()$output,ingest()$standard_info,input$compound_option)}})
        output$rt_table <- renderDT(ingest()$rt_table,options=list('lengthMenu'=JS('[[10,25,50,-1],[10,25,50,"All"]]'),searching=FALSE),class='white-space:nowrap')
        output$data <- renderDT(data()$data,options=list('lengthMenu'=JS('[[10,25,50,-1],[10,25,50,"All"]]'),searching=FALSE),class='white-space:nowrap')
        output$plot_irms_drift <- renderPlot(data()$plot_irms_drift)
        output$plot_area_curve <- renderPlot(data()$plot_area_curve)
        output$plot_ampl_curve <- renderPlot(data()$plot_ampl_curve)
   }
   #Processing Order Tab
    {
        first_correction_data <- reactive({
            if(input$first_correction == "Drift")
                {output = if(is.null(ingest()$output)){NULL}
                          else{drift_function(data()$data,input$drift_option,input$drift_comp)}}
            if(input$first_correction == "Size")
                {output = if(is.null(ingest()$output)){NULL}
                          else{size_function(data()$data,input$size_option,input$size_cutoff,input$size_normal_peak_action,input$size_small_peak_action,
                                             input$size_toosmall_peak_action,
                                             input$size_large_peak_action,input$acceptable_peak_units,input$largest_acceptable_peak,input$smallest_acceptable_peak)}}
            if(input$first_correction == "Scale Normalization")
                {output = if(is.null(ingest()$output)){NULL}
                          else{normalization_function(data()$data,input$normalization_option, input$normalization_comps)}}
            output
        })
        second_correction_data <- reactive({
            if(input$second_correction == "Drift")
                {output = if(is.null(first_correction_data()$output)){NULL}
                          else{drift_function(first_correction_data()$output,input$drift_option,input$drift_comp)}}
            if(input$second_correction == "Size")
                {output = if(is.null(first_correction_data()$output)){NULL}
                          else{size_function(first_correction_data()$output,input$size_option,input$size_cutoff,input$size_normal_peak_action,input$size_small_peak_action,
                                             input$size_toosmall_peak_action,
                                             input$size_large_peak_action,input$acceptable_peak_units,input$largest_acceptable_peak,input$smallest_acceptable_peak)}}
            if(input$second_correction == "Scale Normalization")
                {output = if(is.null(first_correction_data()$output)){NULL}
                          else{normalization_function(first_correction_data()$output,input$normalization_option,input$normalization_comps)}}
            output
        })
        third_correction_data <- reactive({
            if(input$third_correction == "Drift")
                {output = if(is.null(second_correction_data()$output)){NULL}
                          else{drift_function(second_correction_data()$output,input$drift_option,input$drift_comp)}}
            if(input$third_correction == "Size")
                {output = if(is.null(second_correction_data()$output)){NULL}
                          else{size_function(second_correction_data()$output,input$size_option,input$size_cutoff,input$size_normal_peak_action,input$size_small_peak_action,
                                             input$size_toosmall_peak_action,
                                             input$size_large_peak_action,input$acceptable_peak_units,input$largest_acceptable_peak,input$smallest_acceptable_peak)}}
            if(input$third_correction == "Scale Normalization")
                {output = if(is.null(second_correction_data()$output)){NULL}
                          else{normalization_function(second_correction_data()$output,input$normalization_option,input$normalization_comps)}}
            output
        })

        drift_correction <- reactive({
            if(input$first_correction == "Drift"){output = first_correction_data()}
            if(input$second_correction == "Drift"){output = second_correction_data()}
            if(input$third_correction == "Drift"){output = third_correction_data()}
            output
        })

        size_correction <- reactive({
            if(input$first_correction == "Size"){output = first_correction_data()}
            if(input$second_correction == "Size"){output = second_correction_data()}
            if(input$third_correction == "Size"){output = third_correction_data()}
            output
        })

        normalization_correction <- reactive({
            if(input$first_correction == "Scale Normalization"){output = first_correction_data()}
            if(input$second_correction == "Scale Normalization"){output = second_correction_data()}
            if(input$third_correction == "Scale Normalization"){output = third_correction_data()}
            output
        })

        processing_order <- reactive({c(input$first_correction,input$second_correction,input$third_correction)})


        observe({
            second_choice = data.frame(choice = c("Drift","Size","Scale Normalization")) %>%
                filter(choice != input$first_correction) %>%
                .$choice
            updateSelectInput(session,"second_correction",choices = second_choice)
        })

        observe({
            third_choice = data.frame(choice = c("Drift","Size","Scale Normalization")) %>%
                filter(choice != input$first_correction & choice != input$second_correction) %>%
                .$choice
            updateSelectInput(session,"third_correction",choices = third_choice)
        })
    }
   #Drift Correction Tab
    {
        drift_comp_list <- reactive({
            if(!is.null(ingest()$output)){
                data()$data %>% filter(comp != "Ref") %>%  # Removes the reference peaks.
                    filter(grepl("drift",id2)) %>%   # drift identification may be in id2
                    ungroup() %>%
                    unite(comp_class,c(comp,class),sep=" ") %>% select(comp_class) %>% distinct() %>%
                    .$comp_class} else return(NULL)
            })

        observe({
            x = c("Mean drift of all compounds",drift_comp_list())
            updateRadioButtons(session,"drift_comp",choices = x,selected = ingest()$initials$drift_comp)
        })

        output$drift_raw_plot <- renderPlot(drift_correction()$drift_raw_plot)
        output$drift_corrected_standards_plot <- renderPlot(drift_correction()$drift_corrected_standards_plot)
        output$drift_correction_by_row_plot <- renderPlot(drift_correction()$drift_correction_by_row_plot)
        output$drift_corrected_data <- renderDT(drift_correction()$output,options=list('lengthMenu'=JS('[[10,25,50,-1],[10,25,50,"All"]]'),searching=FALSE),class='white-space:nowrap')
    }
   #Size Correction Tab
    {

       output$size_raw_bycomp_plot <- renderPlot(size_correction()$size_raw_bycomp_plot)
       output$size_raw_grouped_plot <- renderPlot(size_correction()$size_raw_grouped_plot)
       output$size_correction_plot <- renderPlot(size_correction()$size_correction_plot)
       output$size_corrected_bycomp_plot <- renderPlot(size_correction()$size_corrected_bycomp_plot)
       output$size_corrected_grouped_plot <- renderPlot(size_correction()$size_corrected_grouped_plot)
       output$size_corrected_data <- renderDT(size_correction()$output,
                                               options = list('lengthMenu' = JS('[[10, 25, 50, -1], [10, 25, 50, "All"]]'),
                                                              searching= FALSE),
                                               class = 'white-space: nowrap',
                                               filter = "top")
   }
   #Scale Normalization Tab
    {
       output$normalization_plot <- renderPlot(normalization_correction()$normalization_plot)
       output$scale_normalization_corrected_data <- renderDT(normalization_correction()$output,
                                                             options = list('lengthMenu' = JS('[[10, 25, 50, -1], [10, 25, 50, "All"]]'),
                                                                            searching= FALSE),
                                                             class = 'white-space: nowrap',
                                                             filter = "top")
   }
   #QA Performance Tab
    {
       control_standards <- reactive({if(!is.null(third_correction_data())){control_function(third_correction_data()$output,input$normalization_comps,processing_order())}
           else return(NULL)})

       output$control_standards <- renderDT(control_standards()$control_standards,
                                                             options = list(dom = 't'),
                                                             class = 'white-space: nowrap')
       output$corrections_plot <- renderPlot(control_standards()$corrections_plot)


   }
   #Sample Output Tab (with Derivatization Correction)
    {

        
        derivatization_results <- reactive({
            derivatization_correction(third_correction_data()$output,
                                      ingest()$derivatization_table,
                                      input$derivatization_option)
        })
        
        output$derivatization_table <- renderDT(ingest()$derivatization_table[,1:4],options=list('lengthMenu'=JS('[[10,25,50,-1],[10,25,50,"All"]]'),searching=FALSE),class='white-space:nowrap')
        output$derivative_dD_selected <- renderTable(derivatization_results()$derivative_dD_selected)
        output$sample_derivatization_table <- renderDT(derivatization_results()$output,options=list('lengthMenu'=JS('[[10,25,50,-1],[10,25,50,"All"]]'),searching=FALSE),class='white-space:nowrap')
        
    }
    
}

shinyApp(ui, server = server)
