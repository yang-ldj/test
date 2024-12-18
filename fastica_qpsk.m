clear;clc;
close all;


load source_data_qpsk.mat 

evm = comm.EVM();

%% 使用基于负熵的FastICA算法的分离  标记为 f1

for i = 1:length(dB_noise)
    %%%%%%%%%% 算法迭代部分 %%%%%%%%%%
    [B,Q,iteration_num_f1] = fastica_achieve1(mixed_data_noised(:,:,i));
    B_FORLOOK(:,:,i) = B;
    Q_FORLOOK(:,:,i) = Q;
    aa(:,:,i) = B'*Q; 
   
    qpsk_ICAedS(:,:,i) = aa(:,:,i)*mixed_data_noised(:,:,i);   %得到分离后的输出信号
    qpsk_rmsEVM_separated_f1(i,:) = evm(source_data_unnoised.',qpsk_ICAedS(:,:,i).');
    qpsk_rmsEVM_separated_f1(i,:) = 20*log10(0.01*qpsk_rmsEVM_separated_f1(i,:));
    rmsEVM_subtract_f1(i,:) = qpsk_rmsEVM_separated_f1(i,:) - qpsk_rmsEVM_unseparated(i,:);
end



%% %%%%%%%% 绘制第1张 alpha-EVM图 %%%%%%%%%%
figure(1);
plot(dB_noise',qpsk_rmsEVM_unseparated(:,2),'-or');hold on;
grid on;
grid minor;
plot(dB_noise',qpsk_rmsEVM_separated_f1(:,2),'-ok');
% plot(dB_noise',rmsEVM_subtract_f1(:,1),'-*m',dB_noise',rmsEVM_subtract_f1(:,2),'-*k');

xlabel('SNR (dB)');ylabel('EVM (dB)');
legend('unseparation gn','f1 separated gn');
% title('xpd = 10dB, alpha = sqrt(0.1),different alpha && EVM');



% % figure(2); scatter(real(mixed_data_noised.'),imag(mixed_data_noised.'));
% figure(1); subplot(121);plot(mixed_data_noised(1,:,12),'.'); hold on;
% plot(source_data_unnoised(1,:),'xr') % axis([-2 2 -1 1]); 
% title('un-ICAed consequence');
% 
% % scatter(real(ICAedS.'),imag(ICAedS.'));
% subplot(122);plot(ICAedS(1,:,16),'.'); hold on;
% plot(source_data_unnoised(1,:),'xr') % axis([-2 2 -1 1]); 
% title('ICAed consequence');
% 
% scatter(real(sn_unnoise),imag(sn_unnoise));


%% %%%%%%%% 绘制第2张 ICA分离前后alpah变化图 %%%%%%%%%%
num = length(dB_noise);
qpsk_aa_2 = zeros(num,2);
qpsk_aa_3 = zeros(2,2,num);
for i = 1:num
%     aa_3(:,:,i) = real (aa(:,:,i)*mix_matrix); 
    qpsk_aa_3(:,:,i) = aa(:,:,i)*mix_matrix;  
    qpsk_aa_2(i,1) = qpsk_aa_3(2,2,i); 
    qpsk_aa_2(i,2) = qpsk_aa_3(2,1,i);

end


figure(2);

% plot(dB_noise',-20*log10(abs(aa_2(:,1))),'-pg',dB_noise',-20*log10(abs(aa_2(:,2))),'-dg');hold on;
% plot(dB_noise',-20*log10(abs(W_f2_2(:,1))),'-pb',dB_noise',-20*log10(abs(W_f2_2(:,2))),'-db');
% plot(dB_noise',-20*log10(abs(W_f3_2(:,1))),'-pm',dB_noise',-20*log10(abs(W_f3_2(:,2))),'-dm');
% grid on; grid minor;
% xlabel('SNR (dB)');ylabel('separated XPD (dB)');
% legend('f1 separated 12','f1 separated 21','f2 separated 12','f2 separated 21','f3 separated 12','f3 separated 21');
% % title('xpd = 10dB alpha=sqrt(0.1), Input & Output alpha difference');


plot(dB_noise',20*log10(abs(qpsk_aa_2(:,1))./abs(qpsk_aa_2(:,2))),'-og');
hold on;
alpha_db = ones(length(dB_noise),1)*alpha_mix;
plot(dB_noise',-20*log10(abs(alpha_db)),'-or')
grid on; grid minor;
xlabel('SNR (dB)');ylabel('separated XPD (dB)');
legend('f1 separated','reference');

% title('xpd = 10dB alpha=sqrt(0.1), Input & Output alpha difference');





%% %%%%%%%% 绘制第3张 alpha-SINR图 %%%%%%%%%%
%计算分离的性能
P_signal = sum((source_data_unnoised.*conj(source_data_unnoised)).');%计算信号的功率P_signal的两个元素分别为sn和gn的功率
%计算SINR
% dB_SINR = zeros(length(dB_noise),2);

for i =1:length(dB_noise)
    % 计算未分离的SINR
    interference_noise_unseparation(:,:,i) = mixed_data_noised(:,:,i) - source_data_unnoised;
    P_ganrao_plus_noise_unseparation(i,:) = sum((interference_noise_unseparation(:,:,i) .*conj(interference_noise_unseparation(:,:,i) )).');    
    qpsk_dB_SINR_unseparation(i,:) = 10*log10(P_signal./P_ganrao_plus_noise_unseparation(i,:));  

    % 计算 分离的SINR    f1
    interference_noise_separation_f1(:,:,i) = qpsk_ICAedS(:,:,i) - source_data_unnoised;
    P_ganrao_plus_noise_separation_f1(i,:) = sum((interference_noise_separation_f1(:,:,i) .*conj(interference_noise_separation_f1(:,:,i) )).');    
    qpsk_dB_SINR_separation_f1(i,:) = 10*log10(P_signal./P_ganrao_plus_noise_separation_f1(i,:)); 

end




figure(3);
plot(dB_noise',qpsk_dB_SINR_unseparation(:,2),'-.pr');
hold on;grid on;
grid minor;
plot(dB_noise',qpsk_dB_SINR_separation_f1(:,2),'-pg');

xlabel('SNR (dB)');ylabel('SINR (dB)');
legend('unseparated gn');
% title('SNR=20dB, different SNR && SINR');


% %% %%%%%%%% 绘制第4张 alpha-迭代次数图 %%%%%%%%%%
% %
% 
% figure(4);
% plot(dB_noise,iteration_num_f1,'-or');
% hold on;grid on;
% grid minor;
% 
% xlabel('SNR (dB)');ylabel('Time of iteration');
% legend('f1');
% % title('different alpha && 迭代次数');

%% %性能指数（PI）
for i = 1:length(dB_noise)
%     PI_unseparation(i) = 1/2*abs(mix_matrix(1,2,i))/abs(mix_matrix(1,1,i)) + ...
%         abs(mix_matrix(2,1,i))/abs(mix_matrix(2,2,i)) + ...
%         abs(mix_matrix(2,1,i))/abs(mix_matrix(1,1,i)) + ...
%         abs(mix_matrix(1,2,i))/abs(mix_matrix(2,2,i)) ;
%    PI_unseparation = 1/2*4*ones(1,length(dB_noise))*0.1^0.5;
     PI_unseparation = 1/2*4*ones(1,length(dB_noise))*abs(alpha_mix);
    qpsk_PI_separation_1(i) = 1/2*(abs(qpsk_aa_3(1,2,i))/abs(qpsk_aa_3(1,1,i)) +abs(qpsk_aa_3(2,1,i))/abs(qpsk_aa_3(2,2,i)) + ...
        abs(qpsk_aa_3(2,1,i))/abs(qpsk_aa_3(1,1,i)) + abs(qpsk_aa_3(1,2,i))/abs(qpsk_aa_3(2,2,i)) );

end

%%%%%%%%%% 绘制第5张 SNR-PI(performance index)图 %%%%%%%%%%
figure(5);
plot(dB_noise,PI_unseparation,'-ko',dB_noise,qpsk_PI_separation_1,'-bp');
hold on;
plot(dB_noise,ones(1,length(dB_noise))*0.01,'-r');
grid on;grid minor;
% xlim([10 30]);ylim([0.05 0.4]);
% title('different SNR & PI');
xlabel('SNR (dB)');ylabel('PI');
legend('PI unseparation','f1 PI separation','performance boundry');




% save qpsk_0.1_20db_data.mat qpsk_ICAedS qpsk_rmsEVM_separated_f1 qpsk_rmsEVM_unseparated ...
%     qpsk_aa_2 ...
%     qpsk_dB_SINR_unseparation qpsk_dB_SINR_separation_f1 ...
%     qpsk_PI_separation_1

figure(6); 

% subplot(121);
plot(mixed_data_noised(1,:,3),'.'); hold on;
plot(source_data_unnoised(1,:),'xr') ;
xlabel('同向分量');ylabel('正交分量');
axis([-2 2 -2 2]); 
% title('XPD=10dB,SNR=15dB,un-ICAed consequence');

figure(7); 
% subplot(122);
plot(qpsk_ICAedS(1,:,3),'.'); hold on;
plot(source_data_unnoised(1,:),'xr') % axis([-2 2 -1 1]); 
% title('XPD=10dB,SNR=15dB,ICAed consequence');
xlabel('同向分量');ylabel('正交分量');
axis([-2 2 -2 2]); 




%%
%%%%%%%%%% 绘制第一张 星座图 %%%%%%%%%%
% 
% figure(1);
% sz = 15;  %用于设定星座图上点的大小
% scatter(real(mixed_data_noised(1,:,11)),imag(mixed_data_noised(1,:,11)),sz,'ro','filled');hold on;
% scatter(real(mixed_data_noised(2,:,11)),imag(mixed_data_noised(2,:,11)),sz,'og','filled');
% scatter(real(sn_unnoise),imag(sn_unnoise),'bo','filled');
% scatter(real(gn_unnoise),imag(gn_unnoise),'^b','filled');
% xlabel('In-Phase');ylabel('Quadrature');
% xlim([-1.2 1.2]);ylim([-1.2 1.2]);
% % title('source signal');
% legend('mixed data noised：sn','mixed data noised：gn','original sn','original gn');
% title('SNR=20dB mixed signal with noise');
% 
% %%%%%%%%%% 绘制第二张 分离后的星座图 %%%%%%%%%%
% figure(2);
% sz = 15;
% scatter(real(ICAedS(1,:,11)'),imag(ICAedS(1,:,11)'),sz,'or','filled');hold on
% scatter(real(ICAedS(2,:,11)'),imag(ICAedS(2,:,11)'),sz,'og','filled');
% scatter(real(sn_unnoise),imag(sn_unnoise),'bo','filled');
% scatter(real(gn_unnoise),imag(gn_unnoise),'^b','filled');
% xlabel('In-Phase');ylabel('Quadrature');
% xlim([-1.2 1.2]);ylim([-1.2 1.2]);
% legend('separated sn','separated gn','original sn','original gn');
% title('alpha = 0.15,SNR=20dB    Separated signal');

% %%%%%%%%%% 绘制第三张 SNR-EVM图 %%%%%%%%%%
% figure(3);
% plot(dB_noise',rmsEVM_subtract_f2(:,1),'-p',dB_noise',rmsEVM_subtract_f2(:,2),'-d');
% xlabel('SNR dB');ylabel('EVM改善情况');
% legend('sn','gn');
% title('alpha = 0.15,different SNR && EVM');
% 
% 
% % MSE
% % mse_matrix = zers(2,5000);
% for i = 1:length(dB_noise)
% error_matrix_unseparation(:,:,i) = mixed_data_noised(:,:,i) - source_data_unnoised;
% error_matrix_separation_f1(:,:,i) = source_data_unnoised - ICAedS(:,:,i);
% mse_value_unseparation(i,:) = sum(transpose(error_matrix_unseparation(:,:,i).*conj(error_matrix_unseparation(:,:,i))))/length(bitdata2);
% mse_value_separation(i,:) = sum(transpose(error_matrix_separation_f1(:,:,i).*conj(error_matrix_separation_f1(:,:,i))))/length(bitdata2);
% end
% %%%%%%%%%% 绘制第四张 SNR-mse图 %%%%%%%%%%
% figure(4);
% plot(dB_noise',mse_value_unseparation(:,1),'-.p',dB_noise',mse_value_unseparation(:,2),'-.d');hold on;
% plot(dB_noise',mse_value_separation(:,1),'-p',dB_noise',mse_value_separation(:,2),'-d');hold on;
% xlabel('SNR dB');ylabel('mse');
% legend('unseparated sn','unseparated gn','separated sn','separated gn');
% title('alpha = 0.15, different SNR && mse');
% 
% 
% %计算分离的性能
% P_signal = sum((source_data_unnoised.*conj(source_data_unnoised)).');%计算信号的功率P_signal的两个元素分别为sn和gn的功率
% %计算SINR
% % dB_SINR = zeros(length(dB_noise),2);
% 
% for i =1:length(dB_noise)
%     % 计算未分离的SINR
%     interference_noise_unsepation(:,:,i) = mixed_data_noised(:,:,i) - source_data_unnoised;
%     P_ganrao_plus_noise_unseparation(i,:) = sum((interference_noise_unsepation(:,:,i) .*conj(interference_noise_unsepation(:,:,i) )).');    
%     dB_SINR_unseparation(i,:) = 10*log10(P_signal./P_ganrao_plus_noise_unseparation(i,:));   
%     % 计算 分离的SINR
%     interference_noise_sepation(:,:,i) = ICAedS(:,:,i) - source_data_unnoised;
%     P_ganrao_plus_noise_separation(i,:) = sum((interference_noise_sepation(:,:,i) .*conj(interference_noise_sepation(:,:,i) )).');    
%     dB_SINR_separation(i,:) = 10*log10(P_signal./P_ganrao_plus_noise_separation(i,:));    
% end
% 
% %%%%%%%%%% 绘制第五张 SNR-SINR图 %%%%%%%%%%
% figure(5);
% plot(dB_noise',dB_SINR_unseparation(:,1),'-.p',dB_noise',dB_SINR_unseparation(:,2),'-.d');hold on;
% plot(dB_noise',dB_SINR_separation(:,1),'-p',dB_noise',dB_SINR_separation(:,2),'-d');hold on;
% xlabel('SNR dB');ylabel('SINR');
% legend('unseparated sn','unseparated gn','separated sn','separated gn');
% title('alpha = 0.15, different SNR && SINR');
% 










