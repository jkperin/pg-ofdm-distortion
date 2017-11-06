clear, clc, close all

load resultados20km

FF = resultados(1).FF;
SIRteo10 = (resultados(1).SIRteo + resultados(2).SIRteo)/2;
SIRteo40 = (resultados(3).SIRteo + resultados(4).SIRteo)/2;
SIRteo75 = (resultados(5).SIRteo + resultados(6).SIRteo)/2;

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
box(axes1,'on');
hold(axes1,'all');

plot1 = plot(FF, SIRteo10, 'ob', FF, resultados(1).SIReva, 'b', FF, resultados(2).SIReva, '--b',...
     FF, SIRteo40, 'ok', FF, resultados(3).SIReva, 'k', FF, resultados(4).SIReva, '--k',...
     FF, SIRteo75, 'or', FF, resultados(5).SIReva, 'r', FF, resultados(6).SIReva, '--r');

xlabel('Frequency (GHz)', 'Fontsize', 14)
ylabel('SIR (dB)', 'Fontsize', 14)
legend('Simulation', 'Proposed theory up to 3rd order', 'Proposed theory up to 4th order')

set(plot1([1 4 7]),'MarkerSize',2,'Marker','o','LineWidth',1)
set(plot1([2 5 8]),'LineWidth',2)
set(plot1([3 6 9]),'LineWidth',2)
set(gca, 'FontSize', 14)
grid on

set(plot1([6]), 'visible', 'off')
set(figure1, 'Position', [100 100 640 550])

% Create textbox
annotation(figure1,'textbox',...
    [0.705935534591199 0.166701048951049 0.219093946540877 0.0747863247863249],...
    'String',{'\kappa P_0 = 75 GHz'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none',...
    'Color',[0 0 1]);

% Create ellipse
annotation(figure1,'ellipse',...
    [0.777161556603779 0.247767857142857 0.0320188679245282 0.169888375513375]);

% Create ellipse
annotation(figure1,'ellipse',...
    [0.779046383647804 0.129090909090909 0.0320188679245282 0.0672727272727273],...
    'Color',[0 0 1]);

% Create textbox
annotation(figure1,'textbox',...
    [0.702330974842772 0.405168359418361 0.221106525157229 0.0747863247863249],...
    'String',{'\kappa P_0 = 40 GHz'},...
    'FontWeight','bold',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'LineStyle','none');

dir = pwd;
cd 'D:\Dropbox\PG\publicações\artigo2\figs'
saveas(gcf, 'OFDM_Gen_20km', 'fig')
saveas(gcf, 'OFDM_Gen_20km', 'epsc')
cd(dir);
