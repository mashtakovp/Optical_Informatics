clear all

a = -3;
b = 3;
p = -3;
q = 3;
n = 1000;
m = 1000;
alpha = 1;
h_x = (b - a)/n;
h_ksi = (q - p)/ m;
x = a:h_x:b - h_x/2;

function f = input_signal(x)
  betta = -1/10;
  f = exp(1i * betta * x);
endfunction


function result = hermite(x)
  n = 5;
  result = zeros(size(x));
  for k = 0:1:n/2
    result = result + power(-1, k) / factorial(n - 2*k) / factorial (k) * power(2*x, n - 2*k);
  end
  result = result * factorial(n);
endfunction


function res_core = core(ksi, x_k, alpha)
  res_core = 1i* exp((-(power(ksi * x_k, 2)))) * hermite(alpha * ksi * x_k);
endfunction

for l=1:1:m
  ksi_l(l) = p + l * h_ksi;
end

A = zeros(m,n);
for l=1:1:m
  for k=1:1:n
    A(l,k)= core(ksi_l(l), x(k), alpha);
  endfor
endfor

inp_f = input_signal(x)
res_integral = A * (inp_f)' * h_x;

function draw(name, x, plotfunc)
  figure();
  plot(x, plotfunc);
  xlabel(name);
endfunction

draw("Амплитуда входного сигнала", x, abs(inp_f));
draw("Фаза входного сигнала", x, arg(inp_f));
draw("Амплитуда выходного сигнала", ksi_l, abs(res_integral));
draw("Фаза выходного сигнала", ksi_l, arg(res_integral));

