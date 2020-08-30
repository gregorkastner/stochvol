svsample_cpp <- function(y, draws, burnin, designmatrix, priorspec, thinpara, thinlatent, keeptime, startpara, startlatent, keeptau, quiet, correct_model_misspecification, interweave, myoffset, fast_sv) {
    .Call(`_stochvol_svsample_cpp`, y, draws, burnin, designmatrix, priorspec, thinpara, thinlatent, keeptime, startpara, startlatent, keeptau, quiet, correct_model_misspecification, interweave, myoffset, fast_sv, PACKAGE = "stochvol")
}

svlsample_cpp <- function(y, draws, burnin, designmatrix, priorspec, thinpara, thinlatent, keeptime, startpara, startlatent, keeptau, quiet, correct_model_misspecification, interweave, myoffset, general_sv) {
    .Call(`_stochvol_svlsample_cpp`, y, draws, burnin, designmatrix, priorspec, thinpara, thinlatent, keeptime, startpara, startlatent, keeptau, quiet, correct_model_misspecification, interweave, myoffset, general_sv, PACKAGE = "stochvol")
}

get_omori_constants <- function() {
    .Call(`_stochvol_get_omori_constants`, PACKAGE = "stochvol")
}

# Register entry points for exported C++ functions
methods::setLoadAction(function(ns) {
    .Call('_stochvol_Export_registerCCallable', PACKAGE = 'stochvol')
})
