close all;
R = 5;
m = 3;
n = 256;
N = 2 * n + 1;
h_r = R / n;
r_l = 0:h_r:(R-h_r / 2);
alpha = 2;
p = 3;
M = 1024;
b = (n ** 2) / (4 * R * M);
h_x = (2 * R) / n;

function result = input_func(rl, p, alpha)
    result = arrayfun(@(r) besselj(abs(p), alpha * r), rl);
end


function standart_draw(arr, draw_function, a_name, p_name)
    figure;
    subplot(1, 2, 1);
    plot(arr, abs(draw_function), 'LineWidth', 2);
    #ylim([-0.1 0.5]);
    #xlim([0 5.5]);
    title(a_name);
    grid on;

    subplot(1, 2, 2);
    plot(arr, angle(draw_function), 'LineWidth', 2);
    #ylim([-0.1 3.5]);
    #xlim([0 5.5]);
    title(p_name);
    grid on;
end

res1 = input_func(r_l, p, alpha);
standart_draw(r_l, res1, "Амплитуда одномерной моды Бесселя", "Фаза одномерной моды Бесселя");

function a = restore(f, N)
    m = 3;
    n = N - 1;
    a = zeros(2 * N + 1, 2 * N + 1);
    for j = 1:(2 * N + 1)
        for k = 1:(2 * N + 1)
            alpha = round(sqrt((j - n - 1) ^ 2 + (k - n - 1) ^ 2));
            if (alpha > n)
              continue;
            else
                a(j, k) = f(alpha + 1) * exp(1i * m * atan2(k - n - 1, j - n - 1));
            end
        end
    end
end

function draw_2d(f, R, a_name, p_name)

    figure;
    subplot(1, 2, 1);
    imagesc([-R R], [-R R], abs(f));
    title(a_name);
    colormap hot;
    colorbar;
    axis square;

    subplot(1, 2, 2);
    imagesc([-R R], [-R R], angle(f));
    title(p_name);
    colormap hot;
    colorbar;
    axis square;

end

res2 = restore(input_func(r_l, p, alpha), n);

draw_2d(res2, R, "Амплитуда восстановленного изображения", "Фаза восстановленного изображения")
function hank_res = hankel(f, arr)

    m = 3;
    x = arr;
    h_r = x(2) - x(1);
    n = length(x);
    hank_res = zeros(n, 1) + 0i;
    for i = 1:n
        x_item = x(i);
        hank_res(i) = sum(f .* besselj(m, 2 * pi * x * x_item) .* x * h_r);
    end
    hank_res = hank_res * (2 * pi / 1i ^ m);
end
tic;
res_hank = hankel(input_func(r_l, p, alpha), r_l);
time1 = toc;
fprintf('Преобразование Ханкеля: %f секунд\n', time1);
res3 = restore(res_hank, n);
standart_draw(r_l, res_hank, "Амплитуда ПХ исходной одномерной функции", "Фаза ПХ исходной одномерной функции")
draw_2d(res3, R, "Амплитуда ПХ востановленной функции", "Фаза ПХ востановленной функции")


function once_fourrier = dft(N, M, h_x, f)
    numZeros = round((M - N) / 2);
    paddedF = [zeros(1, numZeros), f, zeros(1, numZeros)];
    shiftedF = fftshift(paddedF);
    transformedF = fft(shiftedF) * h_x;
    shiftedTransformedF = fftshift(transformedF);
    midPoint = length(shiftedTransformedF) / 2;
    once_fourrier = shiftedTransformedF(midPoint - N / 2 + 1 : midPoint + N / 2);
end

function fourrier = fourrier_2d(N, M, h_x, f)
  fourrier = zeros(N, N);
  for row = 1:N
    fourrier(row, :) = dft(N, M, h_x, f(row, :));
  endfor
  for col = 1:N
    fourrier(:, col) = dft(N, M, h_x, fourrier(:, col).').';
  endfor
end


tic;
double_fft = fourrier_2d(N, M, h_x, res2);
time2 = toc;
fprintf('Преобразование Фурье: %f секунд\n', time2);
draw_2d(double_fft, b, "Амплитуда двумерного ПФ", "Фаза двумерного ПФ")



