clear, clc, close all

%% OFDM GbE
addpath('PGresultados')
load resultados_GbE

% for k= 1:6
%     pteo = polyfit(FF, resultados(k).SIRteo, 4);
%     figure, plot(FF, polyval(pteo, FF), 'k', FF, resultados(k).SIRteo, '*r')
% end

% L = 100.00km, P0 = 5.00mW, mIM = 2.00%, Nit = 1000, kappa = 15.00
% L = 100.00km, P0 = 4.00mW, mIM = 2.00%, Nit = 1000, kappa = 10.00
% L = 100.00km, P0 = 1.00mW, mIM = 2.00%, Nit = 1000, kappa = 10.00

maxkP.pteo = polyfit(FF, resultados(1).SIRteo, 4);
medkP.pteo = [-27.631551581235957 2.420165457619582e+02 -8.662474757795826e+02 1.616428819759815e+03 -1.641364060960882e+03 8.326485224602728e+02 -1.212926726707096e+02]; %polyfit(FF, resultados(2).SIRteo, 4);
% minkP.pteo = polyfit(FF, resultados(3).SIRteo, 4);

maxkP.pap1 = [-24.939930516563372 2.166585577263484e+02 -7.730315997368609e+02 1.447912393469733e+03 -1.489628945043358e+03 7.739485688856221e+02 -1.250905222965702e+02]; %polyfit(FF, resultados(1).SIRap1, 4);
medkP.pap1 = [2.928985958659193 -50.554411499686640 2.851750811608591e+02 -7.671490585541574e+02 1.095274974730229e+03 -8.198152923374641e+02 2.879343850968958e+02];%[-27.734936774296198 2.426730385563699e+02 -8.681252705100085e+02 1.620616241544421e+03 -1.649591073955123e+03 8.416274363035903e+02 -1.245391858141802e+02]; %polyfit(FF, resultados(2).SIRap1, 4);
% minkP.pap1 = polyfit(FF, resultados(3).SIRap1, 4);

maxkP.pap2 = [-25.658658948313290 2.231262142060370e+02 -7.969737578331802e+02 1.494580823063237e+03 -1.540171470682389e+03 8.029138460263696e+02 -1.315638931980600e+02]; %polyfit(FF, resultados(1).SIRap2, 4);
medkP.pap2 = [3.821361258787405 -58.457562390925910 3.141089824355341e+02 -8.232182427760174e+02 1.155949229414481e+03 -8.542972647093054e+02 2.971914621368896e+02]; %[-26.971274030374957 2.367291640415094e+02 -8.495882418884362e+02 1.591496808093112e+03 -1.626388104185168e+03 8.344272522076262e+02 -1.234987997200365e+02]; % polyfit(FF, resultados(2).SIRap2, 4);
% minkP.pap2 = polyfit(FF, resultados(3).SIRap2, 4);

maxkP.SIReva3 = resultados(1).SIReva;
maxkP.SIReva4 = [25.6426190268395 25.4483679176673 25.2578579997669 25.0709023040971 24.8874172554610 24.7072349418763 24.5302870893974 24.3564218554355 24.1855837233616 24.0176330313964 23.8525264754957 23.6901361396527 23.5304262770291 23.3732806448687 23.2186712020147 23.0664890746338 22.9167143218706 22.7692458806835 22.6240679510497 22.4810879153501 22.3402948611252 22.2016008127359 22.0650006307757 21.9304118974004 21.7978315946571 21.6671838968982 21.5384689938063 21.4116140681107 21.2866236944177 21.1634292717507 21.0420362189529 20.9223814444368 20.8044725179461 20.6882483416152 20.5737200127645 20.4608264169796 20.3495837339707 20.2399357286292 20.1318984043194 20.0254168259217 19.9205100445545 19.8171260057430 19.7152831934164 19.6149340812162 19.5160982030864 19.4187289317311 19.3228485294247 19.2284130042579 19.1354436017243 19.0439007704035 18.9538064927470 18.8651218938683 18.7778695307740 18.6920110359834 18.6065342739978 18.5214057987637 18.4376826166937 18.3553319983193 18.2743805576376 18.1947962036862 18.1166038945525 18.0397784162263 17.9643408690663 17.8902687786007 17.8175878561089 17.7462741524766 17.6763536158264 17.6078097313958 17.5406642104698 17.4749035964769 17.4105544852165 17.3476020010384 17.2860729571706 17.2259607145575 17.1672875085984 17.1100502849670 17.0542765964953 16.9999621942353 16.9471349410821 16.8957999457790 16.8459802135926 16.7976852488017 16.7509440266828 16.7057653293904 16.6621787061557 16.6202038674551 16.5798653445644 16.5411884773146 16.5042047539371 16.4689396617673 16.4354258176778 16.4037018876662 16.3737955566069 16.3457530299622 16.3196104894120 16.2954158381479 16.2732074729936 16.2530498266988 16.2349769686446 16.2190639395043 16.2053558328839 16.1939321813604 16.1848424183501 16.1781878755128 16.1740153875703 16.1724376703381 16.1735212377196 16.1773836486276];
maxkP.SIRevaCAD3 = [25.5372053490125 25.3398675212988 25.1463136315500 24.9564369357616 24.7700648955568 24.5871084222976 24.4074133322000 24.2309051966545 24.0574452972251 23.8869714892674 23.7193581980879 23.5545536593588 23.3924435581105 23.2329849714390 23.0760733013726 22.9216732123169 22.7696885496502 22.6200905356233 22.4727903995987 22.3277650711440 22.1849322764925 22.0442739458304 21.9057135548434 21.7692374442532 21.6347742060224 21.5023140966644 21.3717902854231 21.2431965284360 21.1164701114758 20.9916079394147 20.8685510192982 20.7472991085076 20.6277965953451 20.5100458396237 20.3939943178294 20.2796467814165 20.1669535421865 20.0559215661880 19.9465037823817 19.8387092236244 19.7324932481100 19.6278668333532 19.5247876054279 19.4232683870760 19.3232689344952 19.2248038366124 19.1278348631853 19.0323783088388 18.9383978595985 18.8459114724469 18.7548846701562 18.6653370448391 18.5772358931322 18.4906024304668 18.4043669758758 18.3185506292228 18.2341609998598 18.1512210442038 18.0697016820860 17.9896274978078 17.9109710820012 17.8337586806098 17.7579645544713 17.6836166604662 17.6106909446490 17.5392171404874 17.4691729106069 17.4005898482323 17.3334473813441 17.2677790651573 17.2035661607048 17.1408443084705 17.0795966910854 17.0198611811205 16.9616229946680 16.9049224098428 16.8497468144827 16.7961390961718 16.7440889830725 16.6936422113519 16.6447910531861 16.5975843734422 16.5520172330941 16.5081419540062 16.4659566792746 16.4255175731269 16.3868262119370 16.3499430557860 16.3148735349549 16.2816829408141 16.2503810623235 16.2210386569165 16.1936704797385 16.1683535104306 16.1451082043021 16.1240186688516 16.1051119505429 16.0884803761165 16.0741586708314 16.0622487048810 16.0527942185500 16.0459082443537 16.0416451911924 16.0401312495231 16.0414335599073 16.0456939534236 16.0529948996783 16.0634969903046];

medkP.SIReva3 = [33.3447380972556 33.1481647777781 32.9557143202223 32.7673042974395 32.5827205775349 32.4018981500797 32.2246421516806 32.0509020696700 31.8804993177926 31.7133955731307 31.5494261141199 31.3885629648055 31.2306533025549 31.0756780141659 30.9234945588941 30.7740914815762 30.6273351815008 30.4832208769803 30.3416227848303 30.2025419873999 30.0658595730067 29.9315818192144 29.7995958830546 29.6699126824116 29.5424247574621 29.4171472046310 29.2939773583490 29.1729341085628 29.0539190754027 28.9369546211428 28.8219462102990 28.7089194098147 28.5977831439551 28.4885659619887 28.3811799114016 28.2756563401645 28.1719101235792 28.0699752581199 27.9697691869055 27.8713284338860 27.7745727807121 27.6795411835472 27.5861555601464 27.4944572266796 27.4043700582701 27.3159376798161 27.2290857661733 27.1438602189837 27.0601883741945 26.9781183962951 26.8975791609916 26.8186210989817 26.7411745203862 26.6652921422429 26.5898628664900 26.5149385464025 26.4414925272531 26.3695805914659 26.2991353681096 26.2302149928327 26.1627532735439 26.0968107700558 26.0323224205946 25.9693512978296 25.9078334355216 25.8478345267426 25.7892916801999 25.7322733372372 25.6767176758959 25.6226960358678 25.5701476750329 25.5191470060784 25.4696343949040 25.4216875293297 25.3752479311806 25.3303967966188 25.2870768736589 25.2453731354640 25.2052296522030 25.1667354834951 25.1298361475437 25.0946251472123 25.0610496098352 25.0292078936530 24.9990489381077 24.9706764341557 24.9440413865969 24.9192533738413 24.8962657801793 24.8751947187763 24.8559963418965 24.8387940563101 24.8235472622897 24.8103875545447 24.7992781750469 24.7903599666919 24.7836007491278 24.7791518786074 24.7769866681042 24.7772685078262 24.7799773486694 24.7852904568716 24.7931958604066 24.8038869524597 24.8173616601335 24.8338322803067 24.8533089645600 24.8760263747920];
medkP.SIReva4 = [33.3228835029798 33.1269654133771 32.9351580154238 32.7473550316492 32.5633691036060 32.3831116549012 32.2064137640692 32.0332016528220 31.8633219451991 31.6967128550024 31.5332347232677 31.3728364186073 31.2153890858102 31.0608512697112 30.9091038612412 30.7601128544107 30.6137680990806 30.4700425610083 30.3288329327421 30.1901188553182 30.0538034641455 29.9198713604574 29.7882318681353 29.6588745140116 29.5317130869873 29.4067421072523 29.2838797881358 29.1631241802541 29.0443979677540 28.9277029507038 28.8129647568698 28.7001892202167 28.5893050874955 28.4803208857462 28.3731687141780 28.2678589409588 28.1643274138831 28.0625879983627 27.9625781972420 27.8643140652891 27.7677357292400 27.6728619319083 27.5796340195685 27.4880739616171 27.3981247502212 27.3098103000349 27.2230757123281 27.1379474439255 27.0543712919884 26.9723768900336 26.8919112454665 26.8130058706225 26.7356088133190 26.6597534023024 26.5843472337361 26.5094244477985 26.4359740791124 26.3640338174376 26.2935543304589 26.2245752846575 26.1570474003736 26.0910145466975 26.0264266976083 25.9633305430004 25.9016782071686 25.8415178448205 25.7828021543910 25.7255839574049 25.6698149630076 25.6155512301731 25.5627466458005 25.5114590961035 25.4616429890825 25.4133615689993 25.3665681111471 25.3213297424368 25.2776020956103 25.2354547193314 25.1948438599787 25.1558454256517 25.1184145266934 25.0826319237200 25.0484554812915 25.0159693198904 24.9851322321499 24.9560361587307 24.9286389834620 24.9030389681189 24.8791974954475 24.8572176735879 24.8370625056188 24.8188451102543 24.8025282275305 24.7882335865806 24.7759287810523 24.7657428090994 24.7576462749595 24.7517816194513 24.7481206446453 24.7468181387861 24.7478532689116 24.7513922288140 24.7574199736334 24.7661218513000 24.7774871730316 24.7917184954570 24.8088189391748 24.8290081917486];
medkP.SIRevaCAD3 = [32.9734209475161 32.7746451107884 32.5799644904461 32.3893534079642 32.2025415831546 32.0195188365720 31.8400361202226 31.6640960095250 31.4914675039257 31.3221637801778 31.1559693049385 30.9929061506941 30.8327721388313 30.6755968757358 30.5211897912461 30.3695869305886 30.2206078756086 30.0742942244998 29.9304744863764 29.7891950921680 29.6502924404959 29.5138172071862 29.3796127948110 29.2477336430264 29.1180293977189 28.9905578681813 28.8651742863005 28.7419395081182 28.6207137802974 28.5015607417887 28.3843451550441 28.2691332274306 28.1557937992321 28.0443954733767 27.9348107819266 27.8271105860114 27.7211707678557 27.6170643400253 27.5146702315907 27.4140635264830 27.3151259305396 27.2179345425883 27.1223736042076 27.0285221938298 26.9362668733436 26.8456886847995 26.7566763177046 26.6693127796482 26.5834887154589 26.4992891171024 26.4166064309202 26.3355276683812 26.2559469404957 26.1779533293574 26.1004006913605 26.0233757945909 25.9478148025805 25.8738094890674 25.8012568439931 25.7302508303836 25.6606897582448 25.5926698801500 25.5260907485797 25.4610510235189 25.3974514361315 25.3353931905519 25.2747781469209 25.2157102107751 25.1580923362182 25.1020313100170 25.0474311608818 24.9944017617591 24.9428482136371 24.8928837088276 24.8444144371445 24.7975571754778 24.7522192403066 24.7085212949438 24.6663718447721 24.6258957835993 24.5870028963012 24.5498227002945 24.5142663838972 24.4804685386134 24.4483419198458 24.4180267112434 24.3894374466475 24.3627205000842 24.3377924530762 24.3148065640640 24.2936818023532 24.2745791201656 24.2574203021668 24.2423749436273 24.2293681810470 24.2185793730889 24.2099376812744 24.2036335589853 24.1996010390796 24.1980432635383 24.1989002041048 24.2023896151441 24.2084587602527 24.2173423488224 24.2289966612811 24.2436762424906 24.2613486049389 24.2822917083953];

minkP.SIReva3 = resultados(3).SIReva;

fig = figure;
fig1 = plot(FF, polyval(maxkP.pteo, FF), '-k', FF, polyval(maxkP.pap1, FF), '-r', FF, polyval(maxkP.pap2,FF), '-m', FF, maxkP.SIReva3, '*b', FF, maxkP.SIRevaCAD3, 'pg');
hold on
fig2 = plot(FF, polyval(medkP.pteo, FF), '-k', FF, polyval(medkP.pap1, FF), '-r', FF, polyval(medkP.pap2,FF), '-m', FF, medkP.SIReva3, '*b', FF, medkP.SIRevaCAD3, 'pg');
% fig3 = plot(FF, polyval(minkP.pteo, FF), '-k', FF, polyval(minkP.pap1, FF), '-r', FF, polyval(minkP.pap2,FF), '-m', FF, minkP.SIReva3, '*b');
leg = legend('Simulação Sem Aproximações', 'Simulação Log Linear', 'Simulação Log Quadrático', 'Modelo Analítico (3ª ordem)', 'Modelo Analítico CAD (3ª Ordem)');
set(gca, 'FontSize', 12)
set(fig, 'Position', [100 100 640 500])
xlabel('Frequência (GHz)', 'FontSize', 14)
ylabel('SIR (dB)', 'FontSize', 14)
axis([1 2 15 35])
set(fig1(1:3), 'LineWidth', 2)
set(fig2(1:3), 'LineWidth', 2)
% set(fig3(1:3), 'LineWidth', 2)
% set(gca, 'xtick', )
% set(fig3(1:3), 'LineWidth', 2)
grid on
set(leg,'Position',[0.429513888888889 0.693666666666667 0.462326388888889 0.216]);
% Create textbox
annotation(fig,'textbox',[0.676401645768027 0.2005490250826 0.236677115987461 0.0850202429149798],...
    'String',{'\kappa P_0 = 75 GHz'},'FontSize',14,'FitBoxToText','off','LineStyle','none');
% Create ellipse
annotation(fig,'ellipse',[0.6588125 0.48 0.0365 0.192]);
% Create textbox
annotation(fig,'textbox',[0.668560736677121 0.402071292289081 0.236677115987461 0.0850202429149798],...
    'String',{'\kappa P_0 = 40 GHz'},'FontSize',14,'FitBoxToText','off', 'LineStyle','none');
% 
cd('../../DissertacaoPG/figuras/matlab')
saveas(fig, 'OFDM_GbE_100km.fig')
cd ..
saveas(fig, 'OFDM_GbE_100km', 'epsc');
cd('../../codigos/Graficos')
% 
%%
% L = 20.00km, P0 = 5.00mW, mIM = 2.00%, Nit = 1000, kappa = 15.00
% L = 20.00km, P0 = 4.00mW, mIM = 2.00%, Nit = 1000, kappa = 10.00
% L = 20.00km, P0 = 1.00mW, mIM = 2.00%, Nit = 1000, kappa = 10.00

maxkP.pteo = polyfit(FF, resultados(4).SIRteo, 4);
medkP.pteo = [83.077891784441960 -7.793619584521703e+02 3.024267269486824e+03 -6.210870780084659e+03 7.127536554341729e+03 -4.356352906793040e+03 1.170985767370937e+03];
% minkP.pteo = polyfit(FF, resultados(6).SIRteo, 4);
             
maxkP.pap1 = [-17.511342520464847 1.457592671157565e+02 -4.884682412819632e+02 8.348110423357023e+02 -7.425868214122299e+02 2.850531193718631e+02 31.245587719441463]; %polyfit(FF, resultados(4).SIRap1, 4);
medkP.pap1 = [-1.038005226182108e+02 9.459035560646478e+02 -3.553126139456570e+03 7.037734498807211e+03 -7.736883932247023e+03 4.447857552426609e+03 -9.795473611327394e+02]; %[82.528659074489040 -7.743389326434833e+02 3.004802285110274e+03 -6.169016684648910e+03 7.071882765704221e+03 -4.307139119474965e+03 1.144356479679341e+03];
% minkP.pap1 = polyfit(FF, resultados(6).SIRap1, 4);
             
maxkP.pap2 = [-17.566526224752824 1.461609808574889e+02 -4.896256959237618e+02 8.364794479139097e+02 -7.437872010011654e+02 2.855195483064640e+02 32.310642089136046]; %polyfit(FF, resultados(4).SIRap2, 4);
medkP.pap2 = [-1.037209260341250e+02 9.446500230183084e+02 -3.546458747046568e+03 7.020740433073005e+03 -7.714083844658817e+03 4.432458851245462e+03 -9.706063981500176e+02]; %[82.611887634387240 -7.748808909815007e+02 3.006029995590004e+03 -6.169941562736692e+03 7.071658391627291e+03 -4.305697421726653e+03 1.143768024833290e+03];
% minkP.pap2 = polyfit(FF, resultados(6).SIRap2, 4);

maxkP.SIReva3 = resultados(4).SIReva;
maxkP.SIRevaCAD3 = [47.7889849888088 47.5655284579712 47.3463774191582 47.1314174461996 46.9204699236878 46.7134382716070 46.5101622663655 46.3105601849351 46.1144873142462 45.9218744173713 45.7325899762616 45.5465753407542 45.3637103083017 45.1839452801780 45.0071698297108 44.8333421584234 44.6623603437443 44.4941893603661 44.3287347312496 44.1659673556589 44.0057993148231 43.8482067274354 43.6931074845862 43.5404823355690 43.3902543464274 43.2424084033309 43.0968722058843 42.9536343618109 42.8126267405336 42.6738413211661 42.5372137441675 42.4027390644747 42.2703563495556 42.1400634807300 42.0118026549462 41.8855743697782 41.7613236939199 41.6390535648794 41.5187116996923 41.4003033288450 41.2837786241801 41.1691449879862 41.0563548796055 40.9454177748505 40.8362882765199 40.7289778562754 40.6234431371936 40.5196975279156 40.4176995675635 40.3174645603222 40.2189528745676 40.1221816850291 40.0271131188478 39.9337662120756 39.8410659858784 39.7490315237421 39.6586652814141 39.5699885070375 39.4829668458722 39.3976234189560 39.3139254938584 39.2318981022484 39.1515101239975 39.0727885552177 38.9957038932076 38.9202851692807 38.8465045191162 38.7743930981531 38.7039247182711 38.6351327677335 38.5679927907318 38.5025405388799 38.4387533649554 38.3766695388510 38.3162683207607 38.2575906813768 38.2006179127146 38.1453939003128 38.0919021217621 38.0401896278607 37.9902422695798 37.9421105557077 37.8957829383072 37.8513137265430 37.8086942484798 37.7679830147749 37.7291745607221 37.6923320695906 37.6574536828860 37.6246078116924 37.5937966848985 37.5650945983857 37.5385084501024 37.5141192018880 37.4919391263448 37.4720567851831 37.4544906853279 37.4393381117877 37.4266248587411 37.4164582951781 37.4088728001242 37.4039874871292 37.4018469304492 37.4025840317599 37.4062555750510 37.4130107882550 37.4229212076750 37.4361555720058];

medkP.SIReva3 = [57.9793216594455 57.7453880072179 57.5156882033446 57.2901337267552 57.0685213107223 56.8507795301176 56.6367243117884 56.4262984106324 56.2193339790020 56.0157856462644 55.8154993994410 55.6184399015318 55.4244650287218 55.2335479949215 55.0455569675817 54.8604725078520 54.6781717508250 54.4986416209255 54.3217671158545 54.1475407140214 53.9758543465514 53.8067053767017 53.6399918825507 53.4757155565145 53.3137799539410 53.1541906333320 52.9968560539498 52.8417852535504 52.6888911020156 52.5381857928089 52.3895861801695 52.2431073426935 52.0986697492495 51.9562911377165 51.8158952698914 51.6775023551917 51.5410391677425 51.4065282335362 51.2738990938696 51.1431764648415 51.0142924405747 50.8872738258608 50.7620550803172 50.6386650183235 50.5170403015697 50.3972116950565 50.2791179205074 50.1627916530345 50.0481735516996 49.9352981784143 49.8241080246101 49.7146395319963 49.6068369358863 49.5007385665414 49.3952498776626 49.2904085525161 49.1872004665327 49.0856662366523 48.9857534732227 48.8875047267349 48.7908691578445 48.6958913066992 48.6025218557842 48.5108074074513 48.4207001522941 48.3322488448447 48.2454071854237 48.1602261894885 48.0766610849664 47.9947652774926 47.9144955580411 47.8359078743792 47.7589606352636 47.6837125079781 47.6101235953386 47.5382554903231 47.4680700904784 47.3996321535722 47.3329055005445 47.2679583310200 47.2047565504842 47.1433721216570 47.0837732338371 47.0260359858677 46.9701310952810 46.9161392312599 46.8640339379656 46.8139009619770 46.7657170375585 46.7195735831048 46.6754509655799 46.6334469751980 46.5935461512090 46.5558534842593 46.5203583453485 46.4871739120998 46.4562951961514 46.4278447442537 46.4018242054644 46.3783669222729 46.3574824163860 46.3393165611210 46.3238882916377 46.3113581404658 46.3017563927398 46.2952608756392 46.2919156790661 46.2919192219500];
medkP.SIRevaCAD3 = [57.7215808743172 57.4868182250174 57.2562592594473 57.0298773647376 56.8074073297851 56.5888378845032 56.3739248408040 56.1626695330956 55.9548456266075 55.7504649127848 55.5493163439196 55.3514204690332 55.1565794280136 54.9648211717728 54.7759592936678 54.5900280575488 54.4068510651642 54.2264680134525 54.0487112984647 53.8736253359063 53.7010502883558 53.5310347082182 53.3634256484163 53.1982753226021 53.0354369222076 52.8749659333825 52.7167210370201 52.5607606741442 52.4069484516488 52.2553455066698 52.1058198804078 51.9584351965069 51.8130634993647 51.6697707309768 51.5284325593134 51.3891171118706 51.2517033446077 51.1162614680787 50.9826734285739 50.8510114437727 50.7211601855209 50.5931938261157 50.4669995269006 50.3426533832858 50.2200448353310 50.0992518891949 49.9801660753083 49.8628673159614 49.7472490637192 49.6333931791508 49.5211948869343 49.4107380242103 49.3019194544533 49.1948250455270 49.0883114067893 48.9824621607452 48.8782161924319 48.7756619940189 48.6746992456422 48.5754186000462 48.4777210402368 48.3816994823913 48.2872561371050 48.1944863043660 48.1032933606664 48.0137751289800 47.9258361038855 47.8395767909650 47.7549027695692 47.6719174100306 47.5905273586220 47.5108390576440 47.4327602192815 47.3564005932230 47.2816689752467 47.2086786898537 47.1373396552053 47.0677690746623 46.9998780517767 46.9337880150484 46.8694113445961 46.8068740894019 46.7460900305028 46.6871902901098 46.6300902140372 46.5749265181089 46.5216163242971 46.4703025426788 46.4209043404576 46.3735715185182 46.3282256290244 46.2850241756916 46.2438915226373 46.2049938285083 46.1682588050433 46.1338623892201 46.1017363125443 46.0720676249826 46.0447929228337 46.0201119666107 45.9979672834012 45.9785732721572 45.9618797418927 45.9481180772893 45.9372470921437 45.9295180426418 45.9249009601704 45.9236705580597];

minkP.SIReva3 = resultados(6).SIReva;

fig = figure;
fig1 = plot(FF, polyval(maxkP.pteo, FF), '-k', FF, polyval(maxkP.pap1, FF), '-r', FF, polyval(maxkP.pap2,FF), '-m', FF, maxkP.SIReva3, '*b', FF, maxkP.SIRevaCAD3, 'pg');
hold on
fig2 = plot(FF, polyval(medkP.pteo, FF), '-k', FF, polyval(medkP.pap1, FF), '-r', FF, polyval(medkP.pap2,FF), '-m', FF, medkP.SIReva3, '*b', FF, medkP.SIRevaCAD3, 'pg');
% fig3 = plot(FF, polyval(minkP.pteo, FF), '-k', FF, polyval(minkP.pap1, FF), '-r', FF, polyval(minkP.pap2,FF), '-m', FF, minkP.SIReva3, '*b');
leg = legend('Simulação Sem Aproximações', 'Simulação Log Linear', 'Simulação Log Quadrático', 'Modelo Analítico (3ª ordem)', 'Modelo Analítico CAD (3ª Ordem)');
set(gca, 'FontSize', 12)
set(fig, 'Position', [100 100 640 500])
xlabel('Frequência (GHz)', 'FontSize', 14)
ylabel('SIR (dB)', 'FontSize', 14)
% axis([1 2 15 35])
set(fig1(1:3), 'LineWidth', 2)
set(fig2(1:3), 'LineWidth', 2)
% set(fig3(1:3), 'LineWidth', 2)
% set(gca, 'xtick', )
% set(fig3(1:3), 'LineWidth', 2)
grid on
set(leg, 'Position',[0.428789239376524 0.692354925775979 0.463775687913619 0.218623481781377]);

% Create textbox
annotation(fig,'textbox', [0.607623236677122 0.526071292289081 0.236677115987461 0.0850202429149798],...
    'String',{'\kappa P_0 = 40 GHz'}, 'FontSize',14, 'FitBoxToText','off', 'LineStyle','none');
% Create textbox
annotation(fig,'textbox', [0.676401645768027 0.2445490250826 0.236677115987461 0.0850202429149798],...
    'String',{'\kappa P_0 = 75 GHz'}, 'FontSize',14, 'FitBoxToText','off', 'LineStyle','none');
% Create ellipse
annotation(fig,'ellipse',[0.586937500000002 0.454 0.0365 0.124]);
% Create ellipse
annotation(fig,'ellipse',...
    [0.660062500000003 0.190283400809717 0.0365000000000001 0.0870445344129555]);
% 
cd('../../DissertacaoPG/figuras/matlab')
saveas(fig, 'OFDM_GbE_20km.fig')
cd ..
saveas(fig, 'OFDM_GbE_20km', 'epsc');
cd('../../codigos/Graficos')
