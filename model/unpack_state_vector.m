function [h_s,z_s,h_c,H_c,z_c,Vm]=unpack_state_vector(dmesh,v)
h_s=v(1:dmesh.tri.n_elements);
z_s=v(dmesh.tri.n_elements+1:2*dmesh.tri.n_elements);

nn=2*dmesh.tri.n_elements;
h_c=v(nn+1:nn+dmesh.tri.n_edges);
H_c=v(nn+dmesh.tri.n_edges+1:nn+2*dmesh.tri.n_edges);
z_c=v(nn+2*dmesh.tri.n_edges+1:nn+3*dmesh.tri.n_edges);
nn=nn+3*dmesh.tri.n_edges;

Vm=v(nn+1:end);