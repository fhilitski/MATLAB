clear all;
clc;

EI = 25; %pN*um^2
L = 16; %um
r = 0.5; %um
save_img = false; %solve buckling process as a movie

[frames, f, fs, bs, bd, fd] = solve_buckling_BV(EI,L,r,save_img);

figure;
h = plot(fs,f);
hold on;
h = plot(bs,f);
xlabel ('Strain');
ylabel ('Force (pN)');

fd_nm = fd .* 1000;
bd_nm = bd .* 1000;
figure;
h = plot(fd_nm,(f));
hold on;
h = plot(bd_nm,(f));
xlabel ('Distance (\mum)');
ylabel ('Force (pN)');
legend({'filament distance', 'bead distance'});

