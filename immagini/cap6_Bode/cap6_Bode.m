%% bodeCost
%divided in 2 graphs, then composed in illustrator/inkscape
asbode(10, 1,[-2 2],[-45 45],[-225 45])
%uncomment to save file (it will overwrite the old one!)
%print("bodeCost1", '-depsc')
asbode(-10, 1,[-2 2],[-45 45],[-225 45])
%uncomment to save file (it will overwrite the old one!)
%print("bodeCost2", '-depsc') 
%% bodeZPOrig
%composed in illustrator/inkscape
asbode(1, [1 0],[-2 2],[-45 45])
%uncomment to save file (it will overwrite the old one!)
%print("bodeZPOrig1", '-depsc')
%% es1
asbode(1, [1 0 0],[-2 2],[-45 45])
%uncomment to save file (it will overwrite the old one!)
%print("es1", '-depsc')
%% es2
asbode([1 0], 1,[-2 2],[-45 45])
%uncomment to save file (it will overwrite the old one!)
%print("es2", '-depsc')
%% es3
asbode([-2 1], 1,[-2 2],[-45 45])
%uncomment to save file (it will overwrite the old one!)
%print("es3", '-depsc')
%% es4
asbode([1], [2 1],[-2 2],[-45 45])
%uncomment to save file (it will overwrite the old one!)
%print("es4", '-depsc')
%% es5-divided
asbode([100 -100], [-10 9 1],[-3 2],[-45 45],[],0,0,1,0)
%uncomment to save file (it will overwrite the old one!) 
%it will save only phase, the module, must be saved manually
%print("es5-fas", '-depsc') 
%% es5
asbode([100 -100], [-10 9 1],[-3 2],[-45 45],[],0,0,0,0)
%uncomment to save file (it will overwrite the old one!) 
%print("es5", '-depsc') 
%% bodeZPR-Fas polo
%composed in illustrator/inkscape
asbode([1], [1 1],[-2 2],[-45 45],[-120 120],0,0,0,3,0)
%uncomment to save file (it will overwrite the old one!) 
%print("bodeZPR-pol", '-depsc')
%% bodeZPR-Fas zero
%composed in illustrator/inkscape
asbode([1 1],[1] ,[-2 2],[-45 45],[-120 120],0,0,0,3,0)
%uncomment to save file (it will overwrite the old one!) 
%print("bodeZPR-zer", '-depsc') %% bodeZPR-Fas zero
%% Binomio
%composed in illustrator/inkscape
amp = [-90 30]
fas = [-200 20 ]
num = [1]
den = [1 2*1 1]
asbode( num, den,[-2 2],amp,fas,0,0,0,2,0)
% print("binomio/bodeBin-asint1", '-dpdf')
clf
den = [1 2*0.5 1]
asbode( num, den,[-2 2],amp,fas,0,0,0,2,0)
% print("binomio/bodeBin-asint2", '-dpdf')
clf
den = [1 2*0.05 1]
asbode( num, den,[-2 2],amp,fas,0,0,0,2,0)
% print("binomio/bodeBin-asint3", '-dpdf')
clf

opts = bodeoptions('cstprefs');
opts.Grid='on'
opts.Xlim=[0.01 100]

opts.Ylim=amp
opts.PhaseVisible = 'off';
opts.MagVisible = 'on';
for csi=[0:0.2:1],
    den=[1 2*csi 1];
    bodeplot(tf(num,den),opts)
    hold on 
end
%print("binomio/bodeBin-amp", '-dpdf') 
clf

opts.Ylim=fas
opts.PhaseVisible = 'on';
opts.MagVisible = 'off';
for csi=[0:0.2:1],
    den=[1 2*csi 1];
    bodeplot(tf(num,den),opts)
    hold on 
end
%print("binomio/bodeBin-fas", '-dpdf') 

%uncomment to save file (it will overwrite the old one!) 
%% esempio base diagBode
syms x; 
num = sym2poly((-10000*x+1))
den = sym2poly((-10*x+1))
asbode( fliplr(num), fliplr(den),[-1 6],[],[],0,0,0,3,0)
print("diagBode", '-depsc') 
