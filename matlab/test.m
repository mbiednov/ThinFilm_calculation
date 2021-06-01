Num = 70
reflection_S = zeros(Num,1)
reflection_P = zeros(Num,1)

angles = [1:1:Num]

% glass, Film (100 nm), glass.
% angle from 1 to Num (deg).
% wavelength 525 nm.
% the last variable, true means P polarization, false for S.

n_index = [1.5194,1.51,1.5194]
k_index = [0., 0., 0.]
thickness = [0,100,0]

for i0=1:Num
    reflection_S(i0,1) = mexFilm(3,n_index, k_index, thickness,i0,525.0,false)
    reflection_P(i0,1) = mexFilm(3,n_index, k_index, thickness,i0,525.0,true)
end

% read 
x_s=csvread('2um_525nm_S_macleod.csv',1,0)
x_p=csvread('2um_525nm_P_macleod.csv',1,0)

plot(angles, reflection_S*100,'b')
hold on
plot(angles, reflection_P*100,'g')
plot(x_s(:,1), x_s(:,2),'b.')
plot(x_p(:,1), x_p(:,2),'g.')
%xlim([8,71])
hold off
legend('mex S','mex P','macleod S','macloed P')
