close all;
clear;
N = pow2(6);
M = pow2(8);
%N = M;
a = 5;
b = (N .^ 2) / (4 * a * M)
h_x = 2 * a / N;
x = -a:h_x:(a - h_x/2);
h_x_ = 2 * b / N;
x_ = -b:h_x_:(b - h_x_/2);
gauss = exp(-(x.^2));
inp_sig = sin(3*pi*(x));
function F = dft(f, M, h_x)
    N = length(f);
    paddingSize = (M - N) / 2;
    zeroPadding = zeros(1, paddingSize);
    f_padded = [zeroPadding, f, zeroPadding];
    shiftAmount = M/2;
    f_shifted = [f_padded(shiftAmount+1:end), f_padded(1:shiftAmount)];
    F = fft(f_shifted) * h_x;
    F = [F(shiftAmount+1:end), F(1:shiftAmount)];
    F = F(paddingSize + 1:end-paddingSize);
end

inp_sig_bpf = dft(inp_sig, M, h_x);
function f = analytic(x_)
    coefficient = -6 * 1i;
    frequency = 10 * pi;
    denominatorConstant = 9 * pi;
    sinPart = sin(frequency * x_);
    denominator = denominatorConstant - 4 * pi * x_.^2;
    f = (coefficient * sinPart) ./ denominator;
end

function f = analytic2d(x, y)

    coefficient = -36;
    frequency = 10 * pi;
    denominatorConstant = 9 * pi;
    cosPart1 = cos(frequency * (x - y));
    cosPart2 = cos(frequency * (x + y));
    denominatorU = denominatorConstant - 4 * pi * x.^2;
    denominatorV = denominatorConstant - 4 * pi * y.^2;
    f = coefficient * (cosPart1 - cosPart2) ./ (denominatorU .* denominatorV);
end

F = dft(gauss, M, h_x);

figure(1);
plot(x, abs(gauss), "b");
title('Амплитуда Гауссова пучка')

figure(2);
plot(x, arg(gauss), "b");
title('Фаза Гауссова пучка')

figure(3);
plot(x_, abs(F), "b");
title('Амплитуда Гаусс после БПФ')

figure(4);
plot(x_, arg(F), "b");
title('Фаза Гаусс после БПФ')

[X, X_] = meshgrid(x, x_);
core = exp(-2 * pi * 1i * X.* X_);
PF = core * gauss.' * h_x;

figure(5);
plot(x_, abs(PF), "b");
title('Амплитуда Гауссова пучка (финитное преобразование Фурье)')

figure(6);
plot(x_, arg(PF), "b");
title('Фаза гауссова пучка (финитное преобразование Фурье)')

figure(7);
plot(x_, abs(F), "b");

% Первый подграфик
%subplot(1, 2, 1); % Разделяем окно на 1 ряд и 2 колонки и выбираем первую колонку
%plot(x_, abs(F), 'r'); % Рисуем первый график

% Второй подграфик
%subplot(1, 2, 2); % Выбираем вторую колонку той же разметки
%plot(x_, abs(PF), 'b'); % Рисуем второй график

%suptitle('Сравнение амплитуд'); % Общий заголовок для всей фигуры

figure(8);
plot(x_, arg(F), "r");
plot(x_, arg(PF), "b");
%subplot(1, 2, 1);
%plot(x_, arg(F), "r");
%subplot(1, 2, 2);
%plot(x_, arg(PF), "b");
%suptitle('Сравнение фаз')

figure(9); % Создаем новую фигуру с номером 9

% Первый подграфик - Амплитуда входного сигнала
subplot(2, 2, 1); % Разделяем окно на 2 ряда и 2 колонки, выбираем первую ячейку
plot(x, abs(inp_sig), 'b'); % Рисуем график амплитуды входного сигнала
title('Амплитуда входного сигнала');

% Второй подграфик - Фаза входного сигнала
subplot(2, 2, 2); % Выбираем вторую ячейку
plot(x, arg(inp_sig), 'b'); % Рисуем график фазы входного сигнала
title('Фаза входного сигнала');

% Третий подграфик - Амплитуда входного сигнала после БПФ
subplot(2, 2, 3); % Выбираем третью ячейку
plot(x_, abs(inp_sig_bpf), 'r'); % Рисуем график амплитуды входного сигнала после БПФ
title('Амплитуда входного сигнала после БПФ');

% Четвертый подграфик - Фаза входного сигнала после БПФ
subplot(2, 2, 4); % Выбираем четвертую ячейку
plot(x_, arg(inp_sig_bpf), 'r'); % Рисуем график фазы входного сигнала после БПФ
title('Фаза входного сигнала после БПФ');

f = analytic(x_);
%задание 8
figure(10);
subplot(1, 2, 1);
plot(x_, abs(inp_sig_bpf), "r");
title('Амплитуда');
subplot(1, 2, 2);
plot(x_, arg(inp_sig_bpf), "b");
title('Фаза');
%задание 9

y = x.';
f = exp(-(x.^2)-(y.^2));

figure(11);
F = imagesc([-a, a], [-a, a], abs(f));
colormap hot;
colorbar;

figure(12);
F = imagesc([-a, a], [-a, a], arg(f));
colormap hot;
colorbar;

[X, Y] = meshgrid(x, y);
tr = exp(-(X.^2)-(Y.^2));
for row = 1:rows(tr)
 tr(row, :) = dft(tr(row, :), M, h_x);
endfor
for col = 1:columns(tr)
 tr(:, col) = dft(tr(:, col).', M, h_x).';
endfor

figure(13);
F = imagesc([-b, b], [-b, b], abs(tr));
colormap hot;
colorbar;

figure(14);
F = imagesc([-b, b], [-b, b], arg(tr));
colormap hot;
colorbar;


y = x.'
target = sin(3*pi*(x + y));

figure(15);
F = imagesc([-a, a], [-a, a], abs(target));
colormap hot;
colorbar;

figure(16);
F = imagesc([-a, a], [-a, a], arg(target));
colormap hot;
colorbar;

[X, Y] = meshgrid(x, y);
inp_2d = sin(3*pi*(X + Y));
for row = 1:rows(inp_2d)
 inp_2d(row, :) = dft(inp_2d(row, :), M, h_x);
endfor
for col = 1:columns(inp_2d)
 inp_2d(:, col) = dft(inp_2d(:, col).', M, h_x).';
endfor

figure(17);
F = imagesc([-b, b], [-b, b], abs(inp_2d));
colormap hot;
colorbar;

figure(18);
F = imagesc([-b, b], [-b, b], arg(inp_2d));
colormap hot;
colorbar;

[X, Y] = meshgrid(x, y);
F = analytic2d(X, Y);

figure(19);
Res = imagesc([-b, b], [-b, b], abs(F));
colormap hot;
colorbar;

figure(20);
Res = imagesc([-b, b], [-b, b], arg(F));
colormap hot;
colorbar;

