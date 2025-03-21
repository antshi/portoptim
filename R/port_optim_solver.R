#' Numerical Portfolio Optimization with Constraints
#'
#' Calculates the weights of a GMV or tangency portfolio strategy with different constraints and/or penalties.
#'
#' @param sigma_mat a pxp matrix, the covariance matrix of asset returns.
#' @param mu_vec a vector of length p, the expected returns.
#' @param rf a double, the assumed risk-free return. Default value is 0.
#' @param gamma an integer, the assumed risk-averse parameter. Default value is 2.
#' @param Aeq a cxp equality constraint matrix, containing c constraints for p regressors.
#' Default value is Aeq=NULL, no equality constraints.
#' @param beq a cx1 equality constraint vector. Default value is beq=NULL, no equality constraints.
#' @param A a cxp inequality constraint matrix, containing c constraints for p regressors.
#' Default value is A=NULL, no inequality constraints.
#' @param b a cx1 inequality constraint vector. Default value is b=NULL, no inequality constraints.
#' @param lambda1 a tuning parameter value for the lasso penalty. Default value is lambda1=0.
#' @param penidx1 a logical px1 vector,
#' indicating which coefficients are to be penalized with the lasso penalty lambda1.
#' Default value is penidx1=NULL and imposes penalty on all p coefficients.
#' @param lambda2 a tuning parameter value for the ridge penalty. Default value is lambda2=0.
#' @param penidx2 a logical px1 vector,
#' indicating which coefficients are to be penalized with the ridge penalty lambda2.
#' Default value is penidx2=NULL and imposes penalty on all p coefficients.
#' @param gross_c a double, the required gross exposure constraint.
#' Default value is NULL (no constraint). Attention! Works only with GMV.
#' @param porttype a character string.
#' Possible values are "Tang" (tangency portfolio) and "GMV" (global minimum variance portfolio).
#' Default value is "GMV".
#' @param zero_tol a double, indicating the zero tolerance for the calculated weights. Default value is 1e-7.
#' @param res_all a logical.
#' If TRUE, the result includes the calculated weights and the duals from the optimization.
#' If FALSE, only the weights. Default value is FALSE.
#'
#' @details The portfolio optimization with constraints minimizes
#' \deqn{0.5||y - X \beta ||^2_2 + \lambda_1||\beta||_1 + \lambda_2||\beta||^2_2,}
#' subject to \eqn{Aeq \beta = beq} and \eqn{A \beta\le b}.
#'
#' @return a numeric vector with the portfolio weights.
#' @return a list with the portfolio weights and the duals from the optimization.
#'
#' @importFrom ROI OP Q_objective L_constraint ROI_solve
#' @import ROI.plugin.qpoases
#' @importFrom CVXR Variable Minimize quad_form Problem solve
#'
#' @examples
#' data(prices_m)
#' rets_m <- calc_rets(prices_m)
#' sigma_mat <- cov(rets_m)
#' mu_vec <- mu_estim_wrapper(rets_m, mu_estim_ml)
#' port_gmv <- port_optim_gmv(sigma_mat)
#' port_gmv_solver <- port_optim_solver(sigma_mat, mu_vec, porttype = "GMV")
#' port_tang <- port_optim_tang(sigma_mat, mu_vec)
#' port_tang_solver <- port_optim_solver(sigma_mat, mu_vec, porttype = "Tang")
#' @export port_optim_solver
#'
port_optim_solver <-
  function(sigma_mat,
           mu_vec = NULL,
           rf = 0,
           gamma = 2,
           Aeq = NULL,
           beq = NULL,
           A = NULL,
           b = NULL,
           lambda1 = 0,
           penidx1 = NULL,
           lambda2 = 0,
           penidx2 = NULL,
           gross_c = NULL,
           porttype = "GMV",
           zero_tol = 1e-7,
           res_all = FALSE) {
    sigma_mat <- as.matrix(sigma_mat)
    p <- dim(sigma_mat)[1]

    if (is.null(gross_c)) {
      Ident <- diag(1, p)

      # default values penidx1 and penidx2 (all entries are penalized)
      if (is.null(penidx1)) {
        penidx1 <- matrix(TRUE, p, 1)
      }
      dim(penidx1) <- c(p, 1)

      if (is.null(penidx2)) {
        penidx2 <- matrix(TRUE, p, 1)
      }
      dim(penidx2) <- c(p, 1)

      if (porttype == "GMV") {
        # default values Aeq and beq (sum constraint)
        if (is.null(Aeq)) {
          Aeq <- matrix(1, 1, p)
          beq <- rep(1, 1)
        }

        # number of equality constraints
        m1 <- dim(Aeq)[1]

        # default values A and b (no constraint)

        if (is.null(A)) {
          A <- matrix(NA, 0, p)
          b <- rep(0, 0)
        }

        # number of inequality constraints
        m2 <- dim(A)[1]

        # gamma value
        if (gamma != 2) {
          gamma <- 2
          cat("Gamma is different than 2 for GMV. Gamma set to 2.")
        }

        # quadratic coefficient
        H <-
          rbind(
            cbind(
              (gamma / 2) * sigma_mat + lambda2 * Ident,
              -(gamma / 2) * sigma_mat - lambda2 * Ident
            ),
            cbind(
              -(gamma / 2) * sigma_mat - lambda2 * Ident,
              (gamma / 2) * sigma_mat + lambda2 * Ident
            )
          )

        f <- lambda1 * rbind(penidx1, penidx1)

        #  all constraints (first m1 constraints are equality constraints)
        Amatrix <- rbind(cbind(Aeq, -Aeq), cbind(A, -A))
        bvector <- c(beq, b)

        # optimizer
        opt_problem <-
          ROI::OP(
            ROI::Q_objective(H, L = t(f)),
            ROI::L_constraint(
              L = Amatrix,
              dir = c(rep("==", m1), rep("<=", m2)),
              rhs = bvector
            )
          )
        opt <- ROI::ROI_solve(opt_problem, solver = "qpoases")
        opt_sol <- opt$message$primal_solution

        # duals
        duals <- as.numeric(opt$message$dual_solution[-(1:(2 * p))])

        # calculation of actual weights (wpos - wneg)
        betahat <-
          matrix(opt_sol[1:p] - opt_sol[(p + 1):length(opt_sol)], p, 1)

        # round the solutions
        betahat[which(abs(betahat) < zero_tol)] <- 0
      } else if (porttype == "Tang") {
        if (is.null(mu_vec)) {
          stop(cat(
            "Please, enter a value for mu_vec, the expected returns vector."
          ))
        }

        mu_vec_excess <- as.matrix(mu_vec - rf)

        # default values Aeq and beq (no constraint)
        if (is.null(Aeq)) {
          Aeq <- matrix(NA, 0, p)
          beq <- rep(0, 0)
        }

        # number of equality constraints
        m1 <- dim(Aeq)[1]

        # default values A and b (no constraint)

        if (is.null(A)) {
          A <- matrix(NA, 0, p)
          b <- rep(0, 0)
        }

        # number of inequality constraints
        m2 <- dim(A)[1]

        # quadratic coefficient
        H <-
          rbind(
            cbind(
              (gamma / 2) * sigma_mat + lambda2 * Ident,
              -(gamma / 2) * sigma_mat - lambda2 * Ident
            ),
            cbind(
              -(gamma / 2) * sigma_mat - lambda2 * Ident,
              (gamma / 2) * sigma_mat + lambda2 * Ident
            )
          )

        f <-
          rbind(-mu_vec_excess, mu_vec_excess) + lambda1 * rbind(penidx1, penidx1)

        #  all constraints (first m1 constraints are equality constraints)
        Amatrix <- rbind(cbind(Aeq, -Aeq), cbind(A, -A))
        bvector <- c(beq, b)

        # optimizer
        opt_problem <-
          ROI::OP(
            ROI::Q_objective(H, L = t(f)),
            ROI::L_constraint(
              L = Amatrix,
              dir = c(rep("==", m1), rep("<=", m2)),
              rhs = bvector
            )
          )
        opt <- ROI::ROI_solve(opt_problem, solver = "qpoases")
        opt_sol <- opt$message$primal_solution

        # duals
        duals <- as.numeric(opt$message$dual_solution[-(1:(2 * p))])

        # calculation of actual weights (wpos - wneg)
        betahat <-
          matrix(opt_sol[1:p] - opt_sol[(p + 1):length(opt_sol)], p, 1)

        # round the solutions
        betahat[which(abs(betahat) < zero_tol)] <- 0

        # normalization
        if (sum(betahat) != 0) {
          betahat <- betahat / sum(betahat)
        } else {
          betahat <- betahat
        }
      }

      if (res_all) {
        return(list(as.numeric(betahat), as.numeric(duals)))
      } else {
        return(as.numeric(betahat))
      }
    } else {
      w_var <- CVXR::Variable(p)
      obj_func <- CVXR::Minimize(CVXR::quad_form(w_var, sigma_mat))
      constr_gross <- sum(abs(w_var)) <= gross_c
      constr_sum <- sum(w_var) == 1
      problem_form <-
        CVXR::Problem(obj_func, constraints = list(constr_gross, constr_sum))
      result <- CVXR::solve(problem_form)
      weights <- result$getValue(w_var)
      return(as.numeric(weights))
    }
  }
