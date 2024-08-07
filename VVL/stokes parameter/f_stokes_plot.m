function[Eux,Euy]=f_stokes_plot(Ue,Ve,We,n)
Ue=imresize(Ue,[n,n]);
Ve=imresize(Ve,[n,n]);
We=imresize(We,[n,n]);
Z=zeros(size(Ue));
Q=quiver3(Z,Ue,Ve,We,2,'linewidth',1.5);
view(16,82); 
axis tight;
mags = reshape(Q.WData, numel(Q.UData), []);
currentColormap = colormap(jet);
[~, ~, ind] = histcounts(mags, size(currentColormap, 1));
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
set(Q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'
set(Q.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');
 end