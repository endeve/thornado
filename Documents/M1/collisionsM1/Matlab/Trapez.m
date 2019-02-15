function [ Integral ] = Trapez( X, Y )

  Integral = 0.0;
  for i = 1 : size( Y ) - 1
    Integral...
      = Integral...
        + 0.5*(Y(i+1)+Y(i))*(X(i+1)-X(i));
  end

end

