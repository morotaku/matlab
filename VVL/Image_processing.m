clear all; close all; clc;
L=225;X_r=0;Y_r=0;
ang=23;

imp1=importdata('input_file/11 1_0001.ascii.csv');
Imp1=imp1./max(max(imp1));
imp2=importdata('input_file/11 2_0001.ascii.csv');
Imp2=imp2./max(max(imp2));
imp3=importdata('input_file/11 3_0001.ascii.csv');
Imp3=imp3./max(max(imp3));
imp4=importdata('input_file/11 4_0001.ascii.csv');
Imp4=imp4./max(max(imp4));
imp5=importdata('input_file/11 5_0001.ascii.csv');
Imp5=imp5./max(max(imp5));

Imp1_r=rotate(Imp1,ang);
Imp2_r=rotate(Imp2,ang);
Imp3_r=rotate(Imp3,ang);
Imp4_r=rotate(Imp4,ang);
Imp5_r=rotate(Imp5,ang);

figure(1)
subplot(1,5,1);imagesc(Imp1_r);axis equal; axis tight;
subplot(1,5,2);imagesc(Imp2_r);axis equal; axis tight;
subplot(1,5,3);imagesc(Imp3_r);axis equal; axis tight;
subplot(1,5,4);imagesc(Imp4_r);axis equal; axis tight;
subplot(1,5,5);imagesc(Imp5_r);axis equal; axis tight;


Image1=image_cut(Imp1_r,L,X_r,Y_r);
Image2=image_cut(Imp2_r,L,X_r,Y_r);
Image3=image_cut(Imp3_r,L,X_r,Y_r);
Image4=image_cut(Imp4_r,L,X_r,Y_r);
Image5=image_cut(Imp5_r,L,X_r,Y_r);

% Image1=m_filter(Image1,5);
% Image2=m_filter(Image2,5);
% Image3=m_filter(Image3,5);
% Image4=m_filter(Image4,5);
% Image5=m_filter(Image5,5);


Min=0.05;Max=1;
Image1(Image1<=Min)=Min;Image1(Image1>=Max)=Max;Image1=Image1-Min;
Image2(Image2<=Min)=Min;Image2(Image2>=Max)=Max;Image2=Image2-Min;
Image3(Image3<=Min)=Min;Image3(Image3>=Max)=Max;Image3=Image3-Min;
Image4(Image4<=Min)=Min;Image4(Image4>=Max)=Max;Image4=Image4-Min;
Image5(Image5<=Min)=Min;Image5(Image5>=Max)=Max;Image5=Image5-Min;

% Image1=g_filter(Image1,5,0.8);
% Image2=g_filter(Image2,5,0.8);
% Image3=g_filter(Image3,5,0.8);
% Image4=g_filter(Image4,5,0.8);
% Image5=g_filter(Image5,5,0.8);
Image1=m_filter(Image1,3);
Image2=m_filter(Image2,3);
Image3=m_filter(Image3,3);
Image4=m_filter(Image4,3);
Image5=m_filter(Image5,3);

Image1=Image1-Min;Image1_n=Image1./Max;
Image2=Image2-Min;Image2_n=Image2./Max;
Image3=Image3-Min;Image3_n=Image3./Max;
Image4=Image4-Min;Image4_n=Image4./Max;
Image5=Image5-Min;Image5_n=Image5./Max;

figure(2)
colormap('hot')
subplot(1,5,1);imagesc(Image1_n);axis equal; axis tight; axis off;
subplot(1,5,2);imagesc(Image2_n);axis equal; axis tight; axis off;
subplot(1,5,3);imagesc(Image3_n);axis equal; axis tight; axis off;
subplot(1,5,4);imagesc(Image4_n);axis equal; axis tight; axis off;
subplot(1,5,5);imagesc(Image5_n);axis equal; axis tight; axis off;

csvwrite("Image/2/total.csv",Image1_n); %csvで書き出し
csvwrite("Image/2/0.csv",Image2_n); %csvで書き出し
csvwrite("Image/2/45.csv",Image3_n); %csvで書き出し
csvwrite("Image/2/90.csv",Image4_n); %csvで書き出し
csvwrite("Image/2/135.csv",Image5_n); %csvで書き出し

% Center of image
function y=image_cut(image,L,X_r,Y_r)
binary_image = image > 0.5;
image_size=size(binary_image);
col=image_size(1);
lin=image_size(2);

X=0;Y=0;
Sum=sum(sum(binary_image));
for i=1:col
    for j=1:lin
        X=X+j.*binary_image(i,j);
        Y=Y+i.*binary_image(i,j);
    end
end
X_g=int16(X/Sum);
Y_g=int16(Y/Sum);
y=image(Y_g-L+Y_r:Y_g+L+Y_r,X_g-L+X_r:X_g+L+X_r);
end

function y=rotate(img,ang)
    % 回転する角度を指定 (例: 45度)
    rotationAngle = ang;
    
    % 画像のサイズを取得
    [rows, cols, ~] = size(img);
    
    % 回転行列の作成
    theta = deg2rad(rotationAngle);
    rotationMatrix = [cos(theta), -sin(theta); sin(theta), cos(theta)];
    
    % 中心座標の計算
    center = [(cols+1)/2, (rows+1)/2];
    
    % 回転行列を用いて画像を回転
    rotatedImage = zeros(size(img));
    for i = 1:rows
        for j = 1:cols
            % 回転前の座標を計算
            originalCoord = rotationMatrix * ([j; i] - center') + center';
            
            % 回転前の座標が画像の範囲内か確認してから値を代入
            if all(originalCoord > 0) && all(originalCoord <= [cols; rows])
                % バイリニア補間を行う場合は、適切な補間関数を実装してください
                
                rotatedImage(i, j, :) = img(ceil(originalCoord(2)), ceil(originalCoord(1)), :);
            end
        end
    end
    y=rotatedImage;
end

function y=g_filter(img,Size,sigma)
    % ガウシアンフィルタのサイズと標準偏差を指定
    filterSize = Size; % フィルタサイズを調整してください
    sig = sigma;    % ガウシアン分布の標準偏差を調整してください
    
    % ガウシアンカーネルの作成
    [X, Y] = meshgrid(-floor(filterSize/2):floor(filterSize/2), -floor(filterSize/2):floor(filterSize/2));
    gaussianKernel = exp(-(X.^2 + Y.^2) / (2 * sig^2));
    gaussianKernel = gaussianKernel / sum(gaussianKernel(:));
    
    % 画像サイズの取得
    [rows, cols] = size(img);
    
    % 処理結果を格納するための配列
    filteredImage = zeros(size(img));
    
    % ガウシアンフィルタの適用
    for i = floor(filterSize/2)+1:rows-floor(filterSize/2)
        for j = floor(filterSize/2)+1:cols-floor(filterSize/2)
            % フィルターサイズの範囲内の画素値を取得
            neighborhood = img((i-floor(filterSize/2)):(i+floor(filterSize/2)), ...
                                    j-floor(filterSize/2):j+floor(filterSize/2));
            
            % ガウシアンカーネルと画素値の畳み込み
            filteredImage(i,j) = sum(neighborhood(:) .* gaussianKernel(:));
        end
    end
    y=filteredImage;
end

function output_image = m_filter(input_image, window_size)
    [rows, cols] = size(input_image);
    pad_size = floor(window_size / 2);

    output_image = zeros(size(input_image));

    for i = 1:rows
        for j = 1:cols
            row_indices = max(1, i-pad_size):min(rows, i+pad_size);
            col_indices = max(1, j-pad_size):min(cols, j+pad_size);
            window = input_image(row_indices, col_indices); % ウィンドウを取得
            window_vector = window(:); % ウィンドウを1次元ベクトルに変換
            sorted_window = sort(window_vector); % ウィンドウ内のピクセル値をソート
            median_index = floor(numel(sorted_window) / 2) + 1; % メディアンのインデックスを計算
            output_image(i, j) = sorted_window(median_index); % メディアン値を出力画像に代入
        end
    end
end
% function dst = m_filter(img, k)
%     [w, h] = size(img);
%     size = floor(k / 2);
% 
%     % ０パディング処理
%     img = zeros(w + 2*size, h + 2*size);
%     img(size+1:size+w, size+1:size+h) = double(img);
%     dst = img;
% 
%     % フィルタリング処理
%     for x = 1:w
%         for y = 1:h
%             dst(x+size, y+size) = median(img(x:x+k-1, y:y+k-1), 'all');
%         end
%     end
% 
%     dst = dst(size+1:size+w, size+1:size+h);
% end