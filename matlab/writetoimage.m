function V = writetoimage(invec,inimage,num_nodes,outimage)

hdr=spm_vol(inimage);
[img, xyz]=spm_read_vols(hdr);
hdr.fname=outimage;

s=size(img);
img2=zeros(s);

for t=1:num_nodes
  img2(img==t)=invec(t);
end

spm_write_vol(hdr,img2);

V=1;
