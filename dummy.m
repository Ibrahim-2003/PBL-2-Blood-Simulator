function dummy

load('Q_left_heart.mat');

t = 1:length(Q_left_heart);
plot(t,Q_left_heart)

Q = Q_left_heart(1000:end);
mean(Q)

end