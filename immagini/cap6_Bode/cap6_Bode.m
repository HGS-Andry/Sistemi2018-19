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
%composed in illustrator/inkscape
asbode([-2 1], 1,[-2 2],[-45 45])
%uncomment to save file (it will overwrite the old one!)
%print("es3", '-depsc')
%% es4
%composed in illustrator/inkscape
asbode(1, [-2 1],[-2 2],[-45 45])
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