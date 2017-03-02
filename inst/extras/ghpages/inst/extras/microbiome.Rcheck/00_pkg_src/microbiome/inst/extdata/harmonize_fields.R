#' @title Harmonize fields
#' @description Harmonize metadata fields names and contents.
#' @param x data.frame
#' @return data.frame with harmonized fields
#' @examples # m <- harmonize_fields(x)
#' @export
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
harmonize_fields <- function (x) {

  if (("bmi" %in% colnames(x)) && !("bmi_group" %in% colnames(x))) {
    if (is.factor(x$bmi) || is.character(x$bmi)) {
      warning("bmi information is a factor; renaming as bmi_group")
      x$bmi_group <- as.factor(tolower(as.character(x$bmi)))
      x$bmi <- NULL
    }
  }

  if ("bmi_group" %in% colnames(x)) {	 	 
    x$bmi_group <- tolower(x$bmi_group)
    levels <- c("underweight", "lean", "overweight", "obese", "severeobese", "morbidobese", "superobese")
    levels <- levels[levels %in% unique(x$bmi_group)]
    x$bmi_group <- factor(x$bmi_group, levels = levels)
  }

  if ("subject" %in% colnames(x)) {
    x$subject <- factor(x$subject)
  }

  if ("project" %in% colnames(x)) {
    x$project <- factor(x$project)
  }

  if ("sample" %in% colnames(x)) {
    x$sample <- as.character(x$sample)
    if (any(duplicated(x$sample))) {
      stop("duplicated sample IDs - fix !")
    }
  }

  if ("nationality" %in% colnames(x)) {
    x$nationality <- factor(x$nationality)
  }

  if ("group" %in% colnames(x)) {
    x$group <- factor(x$group)
  }

  if (("time" %in% colnames(x)) && !is.numeric(x$time)) {
    x$time <- as.numeric(as.character(x$time))
  }

  if ("sex" %in% colnames(x)) {
    x$sex <- tolower(x$sex)
    x$sex <- gsub("^f$", "female", x$sex)    
    x$sex <- gsub("^m$", "male", x$sex)    
    x$sex <- gsub("^o$", "other", x$sex)    
    levels <- c("female", "male", "other")
    levels <- levels[levels %in% unique(x$sex)]
    x$sex <- factor(tolower(as.character(x$sex)), levels = levels)
  }

  x

}



harmonize_fieldnames <- function (x) {

  x <- harmonize_terms(x, "subject", 
       c("subject", "Subject", "subjectID", "SubjectID", "Subject_ID", "subject_ID"))

  x <- harmonize_terms(x, "sample", 
       c("sample", "Sample", "sampleID", "SampleID", "Sample_ID", "sample_ID"))

  x <- harmonize_terms(x, "project", 
       c("project", "Project", "projectID", "ProjectID", "Project_ID", "project_ID"))

  x <- harmonize_terms(x, "group", 
       c("group", "Group", "groupID", "GroupID", "Group_ID", "group_ID"))

  x <- harmonize_terms(x, "nationality", 
       c("nationality", "Nationality", "nationalityID", "NationalityID", "Nationality_ID", "nationality_ID"))

  x <- harmonize_terms(x, "sex", 
       c("gender", "Gender", "genderID", "GenderID", "Gender_ID", "gender_ID", "sex", "Sex"))

  x <- harmonize_terms(x, "platform", 
       c("platform", "Platform", "SeqTech"))

  x <- harmonize_terms(x, "age", 
       c("age", "Age", "AGE"))

  x <- harmonize_terms(x, "age_group", 
       c("age_group", "Age_group", "age_group"))

  x <- harmonize_terms(x, "bmi", 
       c("bmi", "Bmi"))

  x <- harmonize_terms(x, "bmi_group", 
       c("bmi_group", "Bmi_group", "bmi_group", "BMI_group", "BMI_Group"))

  x <- harmonize_terms(x, "time", 
       c("time", "Time"))

  x <- harmonize_terms(x, "timepoint", 
       c("timepoint", "Timepoint"))

  x <- harmonize_terms(x, "diversity", 
       c("diversity", "Diversity"))

  x <- harmonize_terms(x, "richness", 
       c("richness", "Richness"))

  x

}


harmonize_terms <- function (x, field, synonymes) {

  # Check which synonymes we find
  syno <- intersect(synonymes, x)

  if (length(syno) == 1) {

    # Replace the synonyme if it is unique
    x <- gsub(syno, field, x)

  } else if (sum(synonymes %in% x) > 1) {
    warning(paste("Multiple", field, "fields, harmonize manually."))
  }

  x

}

