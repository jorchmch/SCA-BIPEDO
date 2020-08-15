// Search the SUPERBLOCK in Xcos
for i=1:length(scs_m.objs)
    if typeof(scs_m.objs(i))=="Block" & scs_m.objs(i).gui=="SUPER_f" then
        scs_m = scs_m.objs(i).model.rpar;
        break;
    end
end

// Set the equilibrium point, in this case cruise speed of u=1.5 m/s
X=[1.5;0.001;0.001;0.001];
//X=[0.001;0.001;0.001;0.001];
U=[0.001;0.001];

// linearize the model
sys = lincos(scs_m,X,U);

// obtaingin the matrices A,B,C,D
A=sys.A
B=sys.B
C=sys.C
D=sys.D

// save the data
save("PLANTABIPEDO.sod","X","U","sys")


