

function k_deriv_lij(x::Matrix, log_v0::Float64)
  function(log_l2::Float64)
    k_theta_deriv = k_deriv_logl2(x, log_v0, log_l2)
    k_theta_deriv[1, 2]
  end
end

function k_theta_lij(
  y::Vector,
  x::Matrix,
  log_v0::Float64,
  log_v1::Float64)
  function(log_l2::Float)
    k_theta = kernel(x, x, sqrt(exp(log_l2)), exp(log_v0), exp(log_v1))
    k_theta[1, 2]
  end
end

reference_f = k_theta_lij(y, x, theta[1], theta[2])
test_f = k_deriv_lij(x, theta[1])
check_derivative(reference_f, test_f, 10.0)
check_derivative(reference_f, test_f, 2.0)
check_derivative(reference_f, test_f, 1.0)
check_derivative(reference_f, test_f, 0.0)


function k_deriv_v0(x::Matrix, log_l2::Float64)
  function(log_v0::Float64)
    k_theta_deriv = k_deriv_logv0(x, log_v0, log_l2)
    k_theta_deriv[10, 2]
  end
end

function k_theta_v0(
  y::Vector,
  x::Matrix,
  log_l2::Float64,
  log_v1::Float64)
  function(log_v0::Float64)
    k_theta = kernel(x, x, sqrt(exp(log_l2)), exp(log_v0), exp(log_v1))
    k_theta[10, 2]
  end
end

reference_f = k_theta_v0(y, x, theta[3], theta[2])
test_f = k_deriv_v0(x, theta[3])
check_derivative(reference_f, test_f, 10.0)
check_derivative(reference_f, test_f, 2.0)
check_derivative(reference_f, test_f, 1.0)
check_derivative(reference_f, test_f, 0.0)

function k_deriv_v1(n::Int64)
  function(log_v1::Float64)
    k_theta_deriv = exp(log_v1) * eye(n)
    k_theta_deriv[10, 2]
  end
end

function k_theta_v1(
  y::Vector,
  x::Matrix,
  log_l2::Float64,
  log_v0::Float64)
  function(log_v1::Float64)
    k_theta = kernel(x, x, sqrt(exp(log_l2)), exp(log_v0), exp(log_v1))
    k_theta[10, 2]
  end
end

reference_f = k_theta_v1(y, x, theta[3], theta[1])
test_f = k_deriv_v1(length(y))
check_derivative(reference_f, test_f, 10.0)
check_derivative(reference_f, test_f, 2.0)
check_derivative(reference_f, test_f, 1.0)
check_derivative(reference_f, test_f, 0.0)
