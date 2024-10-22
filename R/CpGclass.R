### =========================================================================
### The CpGStructure class.
### -------------------------------------------------------------------------

setClass("CpGStructure",
         representation(
           CG_islet="GRanges",
           CG_site="GRanges",
           genome="DNAStringSet",
           parameter="list",
		   CG_levels="list",
		   prob="list",
           annot="csAnno",
           Txdb = "TxDb",
           range="list"
         )
)


setGeneric("get_core", function(x) standardGeneric("get_core"))
setGeneric("get_flank", function(x) standardGeneric("get_flank"))
setGeneric("get_linker", function(x) standardGeneric("get_linker"))


setMethod("print", "CpGStructure", function(x){

  CG_islet <- length(x@CG_islet)

  CG_site <- length(x@CG_site)

  types <- table(x@CG_islet$CGislet_type)

  flank <- types[1]

  core <- types[2]

  linker <- types[3]

  min_core_site <- x@parameter$min_core_site

  min_linker_site <- x@parameter$min_linker_site

  Nmin <- x@parameter$Nmin

  Dmax <- x@parameter$Dmax

  M <- x@parameter$M

  CG_var <- paste(x@parameter$CG_var,collapse = ", ")

  cat("CpG Structure object with ",CG_islet," CG islets; ",CG_site," CG sites.\n",sep = "")

  cat("CpG Structure types: ",flank," flank; ",core," core; ",linker," linker.\n",sep = "")

  cat("Parameters: Nmin = ",Nmin,"; Dmax = ",Dmax,"; M = ",M,"; min_core_site = ",min_core_site,
      "; min_linker_site = ",min_linker_site,".\n",sep = "")

  cat("CG variables: ",CG_var,".\n",sep = "")

})

setMethod("show", "CpGStructure",
          function(object){
            CG_islet <- length(object@CG_islet)

            CG_site <- length(object@CG_site)

            types <- table(object@CG_islet$CGislet_type)

            flank <- types[1]

            core <- types[2]

            linker <- types[3]

            min_core_site <- object@parameter$min_core_site

            min_linker_site <- object@parameter$min_linker_site

            Nmin <- object@parameter$Nmin

            Dmax <- object@parameter$Dmax

            M <- object@parameter$M

            CG_var <- paste(object@parameter$CG_var,collapse = ", ")

            cat("CpG Structure object with ",CG_islet," CG islets; ",CG_site," CG sites.\n",sep = "")

            cat("CpG Structure types: ",flank," flank; ",core," core; ",linker," linker.\n",sep = "")

            cat("Parameters: Nmin = ",Nmin,"; Dmax = ",Dmax,"; M = ",M,"; min_core_site = ",min_core_site,
                "; min_linker_site = ",min_linker_site,".\n",sep = "")

            cat("CG variables: ",CG_var,".\n",sep = "")
           })

