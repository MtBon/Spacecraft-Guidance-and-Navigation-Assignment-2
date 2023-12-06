function dy = two_body_rhs(~,y,mu)

r = sqrt(y(1)^2 + y(2)^2 + y(3)^2);

dy = [y(4)
    y(5)
    y(6)
    -(mu/r^3)*y(1)
    -(mu/r^3)*y(2)
    -(mu/r^3)*y(3)];

end