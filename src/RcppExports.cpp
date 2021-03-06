// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace Eigen;

// TransferModalUnfoldingsT
MatrixXd TransferModalUnfoldingsT(MatrixXd S, int d1, int d2, VectorXi dim);
RcppExport SEXP _tensorIMMAM_TransferModalUnfoldingsT(SEXP SSEXP, SEXP d1SEXP, SEXP d2SEXP, SEXP dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type S(SSEXP);
    Rcpp::traits::input_parameter< int >::type d1(d1SEXP);
    Rcpp::traits::input_parameter< int >::type d2(d2SEXP);
    Rcpp::traits::input_parameter< VectorXi >::type dim(dimSEXP);
    rcpp_result_gen = Rcpp::wrap(TransferModalUnfoldingsT(S, d1, d2, dim));
    return rcpp_result_gen;
END_RCPP
}
// Estimation
List Estimation(MatrixXd Y, List Z0, List Sinit, List Ainit, List Binit, List Cinit, List optsList);
RcppExport SEXP _tensorIMMAM_Estimation(SEXP YSEXP, SEXP Z0SEXP, SEXP SinitSEXP, SEXP AinitSEXP, SEXP BinitSEXP, SEXP CinitSEXP, SEXP optsListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< List >::type Z0(Z0SEXP);
    Rcpp::traits::input_parameter< List >::type Sinit(SinitSEXP);
    Rcpp::traits::input_parameter< List >::type Ainit(AinitSEXP);
    Rcpp::traits::input_parameter< List >::type Binit(BinitSEXP);
    Rcpp::traits::input_parameter< List >::type Cinit(CinitSEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    rcpp_result_gen = Rcpp::wrap(Estimation(Y, Z0, Sinit, Ainit, Binit, Cinit, optsList));
    return rcpp_result_gen;
END_RCPP
}
// setuplambda
VectorXd setuplambda(MatrixXd Y, List Z0, List Sinit, List Ainit, List Binit, List Cinit, int nx, int ng, int nlam, VectorXd setlam);
RcppExport SEXP _tensorIMMAM_setuplambda(SEXP YSEXP, SEXP Z0SEXP, SEXP SinitSEXP, SEXP AinitSEXP, SEXP BinitSEXP, SEXP CinitSEXP, SEXP nxSEXP, SEXP ngSEXP, SEXP nlamSEXP, SEXP setlamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< List >::type Z0(Z0SEXP);
    Rcpp::traits::input_parameter< List >::type Sinit(SinitSEXP);
    Rcpp::traits::input_parameter< List >::type Ainit(AinitSEXP);
    Rcpp::traits::input_parameter< List >::type Binit(BinitSEXP);
    Rcpp::traits::input_parameter< List >::type Cinit(CinitSEXP);
    Rcpp::traits::input_parameter< int >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< int >::type ng(ngSEXP);
    Rcpp::traits::input_parameter< int >::type nlam(nlamSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type setlam(setlamSEXP);
    rcpp_result_gen = Rcpp::wrap(setuplambda(Y, Z0, Sinit, Ainit, Binit, Cinit, nx, ng, nlam, setlam));
    return rcpp_result_gen;
END_RCPP
}
// EstPenColumn
List EstPenColumn(MatrixXd Y, List Z0, List Sinit, List Ainit, List Binit, List Cinit, VectorXd lambda, List optsList, List optsList_pen);
RcppExport SEXP _tensorIMMAM_EstPenColumn(SEXP YSEXP, SEXP Z0SEXP, SEXP SinitSEXP, SEXP AinitSEXP, SEXP BinitSEXP, SEXP CinitSEXP, SEXP lambdaSEXP, SEXP optsListSEXP, SEXP optsList_penSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< List >::type Z0(Z0SEXP);
    Rcpp::traits::input_parameter< List >::type Sinit(SinitSEXP);
    Rcpp::traits::input_parameter< List >::type Ainit(AinitSEXP);
    Rcpp::traits::input_parameter< List >::type Binit(BinitSEXP);
    Rcpp::traits::input_parameter< List >::type Cinit(CinitSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    Rcpp::traits::input_parameter< List >::type optsList_pen(optsList_penSEXP);
    rcpp_result_gen = Rcpp::wrap(EstPenColumn(Y, Z0, Sinit, Ainit, Binit, Cinit, lambda, optsList, optsList_pen));
    return rcpp_result_gen;
END_RCPP
}
// EstPenColumnCV
List EstPenColumnCV(MatrixXd Y, List Z0, MatrixXd Ytest, List Ztest0, List Sinit, List Ainit, List Binit, List Cinit, VectorXd lambda, List optsList, List optsList_pen);
RcppExport SEXP _tensorIMMAM_EstPenColumnCV(SEXP YSEXP, SEXP Z0SEXP, SEXP YtestSEXP, SEXP Ztest0SEXP, SEXP SinitSEXP, SEXP AinitSEXP, SEXP BinitSEXP, SEXP CinitSEXP, SEXP lambdaSEXP, SEXP optsListSEXP, SEXP optsList_penSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< List >::type Z0(Z0SEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type Ytest(YtestSEXP);
    Rcpp::traits::input_parameter< List >::type Ztest0(Ztest0SEXP);
    Rcpp::traits::input_parameter< List >::type Sinit(SinitSEXP);
    Rcpp::traits::input_parameter< List >::type Ainit(AinitSEXP);
    Rcpp::traits::input_parameter< List >::type Binit(BinitSEXP);
    Rcpp::traits::input_parameter< List >::type Cinit(CinitSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    Rcpp::traits::input_parameter< List >::type optsList_pen(optsList_penSEXP);
    rcpp_result_gen = Rcpp::wrap(EstPenColumnCV(Y, Z0, Ytest, Ztest0, Sinit, Ainit, Binit, Cinit, lambda, optsList, optsList_pen));
    return rcpp_result_gen;
END_RCPP
}
// EstimationT4
List EstimationT4(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D, List optsList);
RcppExport SEXP _tensorIMMAM_EstimationT4(SEXP YSEXP, SEXP ZSEXP, SEXP SSEXP, SEXP ASEXP, SEXP BSEXP, SEXP CSEXP, SEXP DSEXP, SEXP optsListSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type S(SSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type B(BSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type C(CSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type D(DSEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    rcpp_result_gen = Rcpp::wrap(EstimationT4(Y, Z, S, A, B, C, D, optsList));
    return rcpp_result_gen;
END_RCPP
}
// setuplambdaT4
VectorXd setuplambdaT4(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D, int nlam, VectorXd setlam);
RcppExport SEXP _tensorIMMAM_setuplambdaT4(SEXP YSEXP, SEXP ZSEXP, SEXP SSEXP, SEXP ASEXP, SEXP BSEXP, SEXP CSEXP, SEXP DSEXP, SEXP nlamSEXP, SEXP setlamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type S(SSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type B(BSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type C(CSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type D(DSEXP);
    Rcpp::traits::input_parameter< int >::type nlam(nlamSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type setlam(setlamSEXP);
    rcpp_result_gen = Rcpp::wrap(setuplambdaT4(Y, Z, S, A, B, C, D, nlam, setlam));
    return rcpp_result_gen;
END_RCPP
}
// EstPenColumnT4
List EstPenColumnT4(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D, VectorXd lambda, List optsList, List optsList_pen);
RcppExport SEXP _tensorIMMAM_EstPenColumnT4(SEXP YSEXP, SEXP ZSEXP, SEXP SSEXP, SEXP ASEXP, SEXP BSEXP, SEXP CSEXP, SEXP DSEXP, SEXP lambdaSEXP, SEXP optsListSEXP, SEXP optsList_penSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type S(SSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type B(BSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type C(CSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type D(DSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    Rcpp::traits::input_parameter< List >::type optsList_pen(optsList_penSEXP);
    rcpp_result_gen = Rcpp::wrap(EstPenColumnT4(Y, Z, S, A, B, C, D, lambda, optsList, optsList_pen));
    return rcpp_result_gen;
END_RCPP
}
// EstPenColumnT4CV
List EstPenColumnT4CV(MatrixXd Y, MatrixXd Z, MatrixXd Ytest, MatrixXd Ztest, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D, VectorXd lambda, List optsList, List optsList_pen);
RcppExport SEXP _tensorIMMAM_EstPenColumnT4CV(SEXP YSEXP, SEXP ZSEXP, SEXP YtestSEXP, SEXP ZtestSEXP, SEXP SSEXP, SEXP ASEXP, SEXP BSEXP, SEXP CSEXP, SEXP DSEXP, SEXP lambdaSEXP, SEXP optsListSEXP, SEXP optsList_penSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< MatrixXd >::type Y(YSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type Ytest(YtestSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type Ztest(ZtestSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type S(SSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type A(ASEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type B(BSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type C(CSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type D(DSEXP);
    Rcpp::traits::input_parameter< VectorXd >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< List >::type optsList(optsListSEXP);
    Rcpp::traits::input_parameter< List >::type optsList_pen(optsList_penSEXP);
    rcpp_result_gen = Rcpp::wrap(EstPenColumnT4CV(Y, Z, Ytest, Ztest, S, A, B, C, D, lambda, optsList, optsList_pen));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_tensorIMMAM_TransferModalUnfoldingsT", (DL_FUNC) &_tensorIMMAM_TransferModalUnfoldingsT, 4},
    {"_tensorIMMAM_Estimation", (DL_FUNC) &_tensorIMMAM_Estimation, 7},
    {"_tensorIMMAM_setuplambda", (DL_FUNC) &_tensorIMMAM_setuplambda, 10},
    {"_tensorIMMAM_EstPenColumn", (DL_FUNC) &_tensorIMMAM_EstPenColumn, 9},
    {"_tensorIMMAM_EstPenColumnCV", (DL_FUNC) &_tensorIMMAM_EstPenColumnCV, 11},
    {"_tensorIMMAM_EstimationT4", (DL_FUNC) &_tensorIMMAM_EstimationT4, 8},
    {"_tensorIMMAM_setuplambdaT4", (DL_FUNC) &_tensorIMMAM_setuplambdaT4, 9},
    {"_tensorIMMAM_EstPenColumnT4", (DL_FUNC) &_tensorIMMAM_EstPenColumnT4, 10},
    {"_tensorIMMAM_EstPenColumnT4CV", (DL_FUNC) &_tensorIMMAM_EstPenColumnT4CV, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_tensorIMMAM(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
