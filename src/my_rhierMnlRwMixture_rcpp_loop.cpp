#include "bayesm.h"

//FUNCTION SPECIFIC TO MAIN FUNCTION------------------------------------------------------
int my_rmultinomF(vec p) {
  vec r(1,fill::randu);
  vec cp = cumsum(p);
  int counter=0;
  for(int i=0;i<cp.size();i++){
    if(r[0]>cp[i]){
      counter += 1;
    }
  }
  return(counter+1);
  
  // Rcpp::Environment r_stats("package:stats");
  // Function rmultinom = r_stats["rmultinom"];
};


//[[Rcpp::export]]
double my_llmnl_con(vec const& beta, vec const& y, mat X, int count_out,int initial_price_id,mat initial_price_state_dropped,           //same1
                             int s1,mat transition_matrix_median_steps,mat price_transition_states,           //..
                             vec const& vec_price_states_probs,int draws_length,int number_price_simulations,     //..
                             bool flag_markovian,bool flag_know_state,vec const& SignRes = NumericVector::create(0)){
  
//   // Wayne Taylor 7/8/2016
  
//   // Evaluates log-likelihood for the multinomial logit model WITH SIGN CONSTRAINTS
//   // NOTE: this is exported only because it is used in the shell .R function, it will not be available to users
  
//   //Reparameterize betastar to beta to allow for sign restrictions
//   vec betastar = beta;
  
//   //The default SignRes vector is a single element vector containing a zero
//   //any() returns true if any elements of SignRes are non-zero
//   if(any(SignRes)){ 
//     uvec signInd = find(SignRes != 0);
//     betastar.elem(signInd) = SignRes.elem(signInd) % exp(betastar.elem(signInd));  //% performs element-wise multiplication
//   }
  
//   int n = y.size();
//   int j = X.n_rows/n;
//   mat Xbeta = X*betastar;
  
//   vec xby = zeros<vec>(n);
//   vec denom = zeros<vec>(n);
  
//   for(int i = 0; i<n;i++){      
//     for(int p=0;p<j;p++) denom[i]=denom[i]+exp(Xbeta[i*j+p]);
//     xby[i] = Xbeta[i*j+y[i]-1];
//   }
  
//   return(sum(xby - log(denom)));
// }

 //main
  int n=y.size(); // number of choices
  int j=X.n_rows/n; // number of brands
  int nvar = X.n_cols;
  
  vec Vec(j-1);  //in .r file it was called vec
  
  // added to llmnl function 
  if (count_out != 0) { // if the initial state is not known
    
    // New August 2017 -- Markovian prices ### 
    if (flag_markovian) {
      mat X_forward(j,j+1,fill::zeros);
      for(int i=0;i<j-1;i++){X_forward(i,i)=1;}
      mat dif_probs(number_price_simulations,j,fill::zeros);
      for(int j_markov=0;j_markov<number_price_simulations;j_markov++){
        // draw a price vector
        vec price_draw_this_iter_tmp(draws_length,fill::zeros);
        mat transition_matrix_median_steps_t = transition_matrix_median_steps.t();
        price_draw_this_iter_tmp[0] = my_rmultinomF(transition_matrix_median_steps_t.col(initial_price_id-1));
        for (int i=1;i<draws_length;i++){
          price_draw_this_iter_tmp[i] = my_rmultinomF(transition_matrix_median_steps_t.col(price_draw_this_iter_tmp[i-1]-1));   //check later
        }
        vec price_draw_this_iter(price_draw_this_iter_tmp.size());
        for(int i=0;i<draws_length;i++){
          price_draw_this_iter[i]=price_draw_this_iter_tmp[draws_length-i-1];
        }
        // price_draw_this_iter <- price_draw_this_iter[draws_length:1] //no need anymore
        //# initial random state
        ///////return price_draw_this_iter[0];
        vec temp(j-1);
        temp.fill(1.0/(j-1));
        X_forward(my_rmultinomF(temp)-1,nvar-1) = 1;    //nvar and j+1
        X_forward.col(j-1) = price_transition_states( price_draw_this_iter[0]-1, span(0, j-1) ).t();
        //# forward simulate the choices
        vec prob_forward(j);
        for (int i=1;i<draws_length;i++){
          vec exp_xbeta_forward(j);
          exp_xbeta_forward = exp(X_forward*beta);
          prob_forward = exp_xbeta_forward / sum(exp_xbeta_forward);
          X_forward(my_rmultinomF(prob_forward)-1,nvar-1) = 1;
          //error notes:incompatible matrix dimensions: 10x1 and 1x10
          mat X_forward_tmp = price_transition_states(price_draw_this_iter[i]-1,span(0, j-1));
          X_forward.col(j-1) = X_forward_tmp.t();
        }
        dif_probs.row(j_markov) = prob_forward.t();
        
        //return dif_probs;
      }
      vec tmp_x(j);
      for(int i=0;i<j;i++){    //correspond to apply(,mean)
        tmp_x[i]=mean(dif_probs.col(i));
      }
      //# marginal initial state probability vector
      Vec = tmp_x(span(0,j-1-1))/sum(tmp_x(span(0,j-1-1)));
      
    } else {
        //####### Older iid price case #########
      mat xbeta(j,price_transition_states.n_rows);    //tiny change
      vec tmp_beta = beta(span(0,j-1));  tmp_beta[j-1] = 0;
      mat tmpmat = price_transition_states.cols(0,j-1) * beta[j-1];
      tmpmat = tmpmat.t();
      for(int i=0;i<j;i++){xbeta.row(i) = tmp_beta[i]+tmpmat.row(i);}
      xbeta = xbeta.t();
      mat Prob0(j-1,j-1,fill::zeros);
      for(int j_state=0;j_state<j-1;j_state++) {
        mat xbeta_state=xbeta;
        xbeta_state.col(j_state) = xbeta.col(j_state) + beta[nvar-1];   //error! change nvar to nvar-1
        mat xb0(xbeta_state.n_rows,xbeta_state.n_cols);
        xb0 = exp(xbeta_state);
        mat prob0(xb0.n_rows,xb0.n_cols);
        for(int i=0;i<xb0.n_rows;i++){prob0.row(i) = xb0.row(i)/sum(xb0.row(i));}  // because nan exists in xb0
        //prob0: same size xb0 or xbeta_state or xbeta or t(price_transition_states); 
        vec tmp_prob(prob0.n_cols) ;
        tmp_prob = (vec_price_states_probs.t() * prob0).t();    //careful!!
        Prob0.col(j_state) =  tmp_prob(span(0,j-1-1)); //# transposed transition matrix
        Prob0(j_state,j_state) = Prob0(j_state,j_state) + tmp_prob[j-1]; //# if choose outside option, state stays the same
      }
      
      cx_vec eigval;
      cx_mat eigvec;
      eig_gen(eigval, eigvec, Prob0);
      Vec=abs(eigvec.col(0));
      Vec = real(Vec);
      Vec = Vec / sum(Vec);
      //# these are marginal probabilites of chosings states for average prices
    }
    if (!flag_know_state) {// # if we are in the case when the data was not dropped and initial state is not known
      int s0 = my_rmultinomF(Vec); //# draw for the initial state
      
      for(int i=0;i<count_out*j;i++) {
        
        X(i,nvar-1) = 0;     //due to rewrite X, cannot input type const& X
        
      }
      
      for(int i=s0;i<=s0+(count_out-1)*j;i=i+j) {
        X(i-1,nvar-1) = 1;  //# assign drawn initial state to all observations until the first brand is chosen
      }
    }
  }
  
  //######################################################################################################
  mat Xbeta_tmp = X*beta;
  Xbeta_tmp.reshape(j,Xbeta_tmp.n_cols*Xbeta_tmp.n_rows/j);
  mat Xbeta = Xbeta_tmp.t();    //size: n*j;  row=n??
  //Xbeta=matrix(Xbeta,byrow=T,ncol=j)
  
  //ind=cbind(c(1:n),y)
  vec xby(n);
  for(int i=0;i<n;i++){xby[i] = Xbeta(i,y[i]-1);}
  Xbeta=exp(Xbeta);
  vec iota(j,fill::ones);
  vec denom=log(Xbeta*iota); //n*1
  
  if (flag_know_state) {
    mat xb_init_tmp(initial_price_state_dropped.n_rows,1);
    xb_init_tmp = exp(initial_price_state_dropped*beta); //matrix of expxbeta by each possible (m-1) state
    //xb_init = matrix(exp(initial_price_state_dropped%*%beta), ncol = j, byrow = T) //later reshape
    xb_init_tmp.reshape(j,xb_init_tmp.n_cols*xb_init_tmp.n_rows/j);
    mat xb_init = xb_init_tmp.t();
    mat prob_init(xb_init.n_rows,xb_init.n_cols); // matrix of choice probabilities
    for(int i=0;i<xb_init.n_rows;i++){prob_init.row(i) = xb_init.row(i)/sum(xb_init.row(i));}
    mat vec_with_init_tmp(1,prob_init.n_cols);
    vec_with_init_tmp = Vec.t()*prob_init; //marginal probabilities
    vec vec_with_init = vec_with_init_tmp.t();
    return(sum(xby-denom) + log(vec_with_init[s1-1]));
  } else {
    return(sum(xby-denom));
  }
 }

//mnlRwMetropOnce=
//function(y,X,oldbeta,oldll,s,inc.root,betabar,rootpi){ 
//#
//# function to execute rw metropolis for the MNL
//# y is n vector with element = 1,...,j indicating which alt chosen
//# X is nj x k matrix of xvalues for each of j alt on each of n occasions
//# RW increments are N(0,s^2*t(inc.root)%*%inc.root)
//# prior on beta is N(betabar,Sigma)  Sigma^-1=rootpi*t(rootpi)
//#  inc.root, rootpi are upper triangular
//#  this means that we are using the UL decomp of Sigma^-1 for prior 
//# oldbeta is the current
//     stay=0
//     betac=oldbeta + s*t(inc.root)%*%(matrix(rnorm(ncol(X)),ncol=1))
//     cll=llmnl(betac,y,X)
//     clpost=cll+lndMvn(betac,betabar,rootpi)
//     ldiff=clpost-oldll-lndMvn(oldbeta,betabar,rootpi)
//     alpha=min(1,exp(ldiff))
//     if(alpha < 1) {unif=runif(1)} else {unif=0}
//     if (unif <= alpha)
//             {betadraw=betac; oldll=cll}
//           else
//             {betadraw=oldbeta; stay=1}
//return(list(betadraw=betadraw,stay=stay,oldll=oldll))
//}

mnlMetropOnceOut mnlMetropOnce_con_initialcond(vec const& y, mat X, vec const& oldbeta,                    //sameI     new
                                   double oldll,double s, mat const& incroot,                       //..
                                   vec const& betabar, mat const& rootpi,                       //sameIII     new
                                   int count_out,int initial_price_id,mat initial_price_state_dropped,     //same1
                                   int s1,mat const& transition_matrix_median_steps,mat const& price_transition_states, //..
                                   vec const& vec_price_states_probs,int draws_length,int number_price_simulations,   //..
                                   bool flag_markovian,bool flag_know_state,vec const& SignRes = NumericVector::create(2)){ 
  // Wayne Taylor 10/01/2014
  
  // function to execute rw metropolis for the MNL
  // y is n vector with element = 1,...,j indicating which alt chosen
  // X is nj x k matrix of xvalues for each of j alt on each of n occasions
  // RW increments are N(0,s^2*t(inc.root)%*%inc.root)
  // prior on beta is N(betabar,Sigma)  Sigma^-1=rootpi*t(rootpi)
  //  inc.root, rootpi are upper triangular
  //  this means that we are using the UL decomp of Sigma^-1 for prior 
  // oldbeta is the current
  
  
  mnlMetropOnceOut out_struct;
  

  double unif;
  vec betadraw, alphaminv;
  
  int stay = 0;
  vec betac = oldbeta + s*trans(incroot)*as<vec>(rnorm(X.n_cols));
  double cll = my_llmnl_con(betac,y,X,count_out,initial_price_id,initial_price_state_dropped,
                                               s1,transition_matrix_median_steps,price_transition_states,
                                               vec_price_states_probs,draws_length,number_price_simulations,
                                               flag_markovian,flag_know_state,SignRes);
  double clpost = cll+lndMvn(betac,betabar,rootpi);
  double ldiff = clpost-oldll-lndMvn(oldbeta,betabar,rootpi);
  alphaminv << 1 << exp(ldiff);
  double alpha = min(alphaminv);
  
  if(alpha < 1) {
    unif = as_scalar(vec(runif(1)));
  } else { 
    unif=0;}
  if (unif <= alpha) {
    betadraw = betac;
    oldll = cll;
  } else {
    betadraw = oldbeta;
    stay = 1;
  }
  
  out_struct.betadraw = betadraw;
  out_struct.stay = stay;  
  out_struct.oldll = oldll;
  
  return (out_struct);
}

//MAIN FUNCTION-------------------------------------------------------------------------------------

//[[Rcpp::export]]
List my_rhierMnlRwMixture_rcpp(List const& lgtdata, mat const& Z,                           //SameA     new
                                 vec const& deltabar, mat const& Ad, mat const& mubar, mat const& Amu,          //..
                                 double nu, mat const& V, double s,                           //..
                                 int R, int keep, int nprint, bool drawdelta,                     //..
                                 mat olddelta,  vec const& a, vec oldprob, mat oldbetas, vec ind, vec const& SignRes){ //SameE     new
                                 //int count_out,int initial_price_id,
                                 // Rcpp::Nullable<mat> initial_price_state_dropped = R_NilValue,
                                 // Rcpp::Nullable<int> s1 = R_NilValue){
                                 
                                 // mat const& initial_price_state_dropped,       //same1
                                 // int s1){

                                 // mat const& transition_matrix_median_steps,mat const& price_transition_states, //..
                                 // vec const& vec_price_states_probs,int draws_length,int number_price_simulations,   //..
                                 // bool flag_markovian,bool flag_know_state){

// Wayne Taylor 10/01/2014

  int nlgt = lgtdata.size();
  int nvar = V.n_cols;
  int nz = Z.n_cols;
  
  mat rootpi, betabar, ucholinv, incroot;
  int mkeep;
  mnlMetropOnceOut metropout_struct;
  List lgtdatai, nmix;
  List lgtdata1 = lgtdata[0];

  ///////
  mat transition_matrix_median_steps = as<mat>(lgtdata1["transition_matrix_median_steps"]);
  mat price_transition_states = as<mat>(lgtdata1["price_transition_states"]);
  vec vec_price_states_probs = as<vec>(lgtdata1["vec_price_states_probs"]);
  int draws_length = as<int>(lgtdata1["draws_length"]);
  int number_price_simulations = as<int>(lgtdata1["number_price_simulations"]);
  bool flag_markovian = as<bool>(lgtdata1["flag_markovian"]);
  bool flag_know_state = as<bool>(lgtdata1["flag_know_state"]);
  ///////

  // convert List to std::vector of struct
  std::vector<moments> lgtdata_vector;
  moments lgtdatai_struct;
  for (int lgt = 0; lgt<nlgt; lgt++){
    lgtdatai = lgtdata[lgt];
    
    lgtdatai_struct.y = as<vec>(lgtdatai["y"]);
    lgtdatai_struct.X = as<mat>(lgtdatai["X"]);
    lgtdatai_struct.count_out = as<int>(lgtdatai["count_out"]);
    lgtdatai_struct.initial_price_id = as<int>(lgtdatai["initial_price_id"]);
    lgtdatai_struct.hess = as<mat>(lgtdatai["hess"]);
    // Rcpp::Nullable<mat> initial_price_state_dropped = lgtdatai["initial_price_state_dropped"];
    // Rcpp::Nullable<int> s1 = lgtdatai["s1"];
    // if(initial_price_state_dropped.isNotNull()){
    //   lgtdatai_struct.initial_price_state_dropped = as<mat>(lgtdatai["initial_price_state_dropped"]);
    // }else{
    //   lgtdatai_struct.initial_price_state_dropped = mat(1,1,fill::ones);
    // }
    // if(s1.isNotNull()){
    //   lgtdatai_struct.s1 = as<int>(lgtdatai["s1"]);
    // }else{
    //   lgtdatai_struct.s1 = 1;
    // }
    lgtdatai_struct.initial_price_state_dropped = as<mat>(lgtdatai["initial_price_state_dropped"]);
    lgtdatai_struct.s1 = as<int>(lgtdatai["s1"]);
    lgtdata_vector.push_back(lgtdatai_struct);    
  }
    
  // allocate space for draws
  vec oldll = zeros<vec>(nlgt);
  cube betadraw(nlgt, nvar, R/keep);
  mat probdraw(R/keep, oldprob.size());
  vec loglike(R/keep);
  mat Deltadraw(1,1); if(drawdelta) Deltadraw.zeros(R/keep, nz*nvar);//enlarge Deltadraw only if the space is required
  List compdraw(R/keep);
  
  if (nprint>0) startMcmcTimer();
    
  for (int rep = 0; rep<R; rep++){
    
    //first draw comps,ind,p | {beta_i}, delta
    // ind,p need initialization comps is drawn first in sub-Gibbs
    List mgout;
    if(drawdelta) {
      olddelta.reshape(nvar,nz);
      mgout = rmixGibbs (oldbetas-Z*trans(olddelta),mubar,Amu,nu,V,a,oldprob,ind);
    } else {
      mgout = rmixGibbs(oldbetas,mubar,Amu,nu,V,a,oldprob,ind);
    }
    
    List oldcomp = mgout["comps"];
    oldprob = as<vec>(mgout["p"]); //conversion from Rcpp to Armadillo requires explict declaration of variable type using as<>
    ind = as<vec>(mgout["z"]);
    
    //now draw delta | {beta_i}, ind, comps
    if(drawdelta) olddelta = drawDelta(Z,oldbetas,ind,oldcomp,deltabar,Ad);
    
    //loop over all LGT equations drawing beta_i | ind[i],z[i,],mu[ind[i]],rooti[ind[i]]
      for(int lgt = 0; lgt<nlgt; lgt++){
        List oldcomplgt = oldcomp[ind[lgt]-1];
        rootpi = as<mat>(oldcomplgt[1]);
        
        //note: beta_i = Delta*z_i + u_i  Delta is nvar x nz
        if(drawdelta){
          olddelta.reshape(nvar,nz);
          betabar = as<vec>(oldcomplgt[0])+olddelta*vectorise(Z(lgt,span::all));
        } else {
          betabar = as<vec>(oldcomplgt[0]);
        }
        
        if (rep == 0) oldll[lgt] = my_llmnl_con(vectorise(oldbetas(lgt,span::all)),lgtdata_vector[lgt].y,lgtdata_vector[lgt].X,
                                                        lgtdata_vector[lgt].count_out,lgtdata_vector[lgt].initial_price_id,lgtdata_vector[lgt].initial_price_state_dropped,
                                                        lgtdata_vector[lgt].s1,transition_matrix_median_steps,price_transition_states,
                                                        vec_price_states_probs,draws_length,number_price_simulations,
                                                        flag_markovian,flag_know_state,SignRes);
        
        //compute inc.root
        ucholinv = solve(trimatu(chol(lgtdata_vector[lgt].hess+rootpi*trans(rootpi))), eye(nvar,nvar)); //trimatu interprets the matrix as upper triangular and makes solve more efficient
        incroot = chol(ucholinv*trans(ucholinv));
                
        metropout_struct = mnlMetropOnce_con_initialcond(lgtdata_vector[lgt].y,lgtdata_vector[lgt].X,vectorise(oldbetas(lgt,span::all)),
                                           oldll[lgt],s,incroot,betabar,rootpi,
                                           lgtdata_vector[lgt].count_out,lgtdata_vector[lgt].initial_price_id,lgtdata_vector[lgt].initial_price_state_dropped,
                                           lgtdata_vector[lgt].s1,transition_matrix_median_steps,price_transition_states,
                                           vec_price_states_probs,draws_length,number_price_simulations,
                                           flag_markovian,flag_know_state,SignRes);
         
         oldbetas(lgt,span::all) = trans(metropout_struct.betadraw);
         oldll[lgt] = metropout_struct.oldll;  
      }
      
    //print time to completion and draw # every nprint'th draw
    if (nprint>0) if ((rep+1)%nprint==0) infoMcmcTimer(rep, R);
    
    if((rep+1)%keep==0){
      mkeep = (rep+1)/keep;
      betadraw.slice(mkeep-1) = oldbetas;
      probdraw(mkeep-1, span::all) = trans(oldprob);
      loglike[mkeep-1] = sum(oldll);
      if(drawdelta) Deltadraw(mkeep-1, span::all) = trans(vectorise(olddelta));
      compdraw[mkeep-1] = oldcomp;
    }
  }
  
  if (nprint>0) endMcmcTimer();
  
  nmix = List::create(Named("probdraw") = probdraw,
    		  Named("zdraw") = R_NilValue, //sets the value to NULL in R
				  Named("compdraw") = compdraw);
  
  //ADDED FOR CONSTRAINTS
  //If there are sign constraints, return f(betadraws) as "betadraws"
  //conStatus will be set to true if SignRes has any non-zero elements
  bool conStatus = any(SignRes);
  
  if(conStatus){
    int SignResSize = SignRes.size();
    
    //loop through each sign constraint
    for(int i = 0;i < SignResSize; i++){
      
      //if there is a constraint loop through each slice of betadraw
      if(SignRes[i] != 0){
        for(int s = 0;s < (R/keep); s++){
          betadraw(span(),span(i),span(s)) = SignRes[i] * exp(betadraw(span(),span(i),span(s)));
        }
      }
      
    }//end loop through SignRes
  }
  
  if(drawdelta){
    return(List::create(
        Named("Deltadraw") = Deltadraw,
        Named("betadraw") = betadraw,
        Named("nmix") = nmix,
        Named("loglike") = loglike,
        Named("SignRes") = SignRes));  
  } else {
    return(List::create(
        Named("betadraw") = betadraw,
        Named("nmix") = nmix,
        Named("loglike") = loglike,
        Named("SignRes") = SignRes));
  }
  
}
