function draw3sigma(Ps, ms)
  t = linspace(0, 2*pi, 40);
  N = 100;
  for a=1:size(Ps, 2)
    S = chol(Ps{a}, 'lower'); %dekompozice: SS' = P
    x_num = [];
    t = linspace(0, 2*pi, N);
    for i=1:N
      x_num(:,i) = 3*S*[cos(t(i)); sin(t(i))] + ms{a};  %x = 3*S*u + b(a);
    end
    plot(x_num(1,:), x_num(2,:), 'r-')
  end
end