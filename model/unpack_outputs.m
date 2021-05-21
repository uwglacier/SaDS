function model_output=unpack_outputs(dmesh,yout)
% Function unpack_outputs unpacks model outputs from the matlab
% time-stepper into a friendly structure with fields corresponding to model
% state variables.
% model_output = 
% 
%   struct with fields:
% 
%     hs: [n_elements x n_outputs double]
%     zs: [n_elements x n_outputs double]
%     hc: [n_edges x n_outputs double]
%     Hc: [n_edges x n_outputs double]
%     zc: [n_edges x n_outputs double]
%     Vm: [n_moulins x n_outputs double]


nouts=size(yout,1);
model_output.hs=zeros(dmesh.tri.n_elements,nouts);
model_output.zs=zeros(dmesh.tri.n_elements,nouts);
model_output.hc=zeros(dmesh.tri.n_edges,nouts);
model_output.Hc=zeros(dmesh.tri.n_edges,nouts);
model_output.zc=zeros(dmesh.tri.n_edges,nouts);
% model_output.Vm=zeros(dmesh.tri.n_n

for ii=1:nouts
   [hs,zs,hc,Hc,zc,vm]=unpack_state_vector(dmesh,yout(ii,:)');
   model_output.hs(:,ii)=hs;
   model_output.zs(:,ii)=zs;
   model_output.hc(:,ii)=hc;
   model_output.Hc(:,ii)=Hc;
   model_output.zc(:,ii)=zc;
   model_output.Vm(:,ii)=vm;
end