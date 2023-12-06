function residual = costfunction(x, t_span, W_m, meas_real, mu)

residual = zeros(size(meas_real)); % Initialize output variable

% Propagate x to the epochs of the measurements
fun = @(t,x) two_body_rhs(t,x,mu);
options = odeset('Reltol',1.e-13,'Abstol',1.e-20);
[~,x_prop] = ode113(fun,t_span,x,options);

% Compute predicted measurements (in this case we simulate to have
% directly measured the state)
meas_pred = x_prop;

% Compute the residual of the measurements and append it to the output
for k=1:length(t_span)
    diff_meas_weighted = W_m * (meas_pred(k,:) - meas_real(k,:))';
    residual(k,:) = diff_meas_weighted';
end

end
