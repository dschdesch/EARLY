function [] = plot_monaural_coincidence

CIRevdir = 'C:\work\BinRev\CIRev';
cd (CIRevdir);
filenames_across = dir('*_across_*');
filenames_within = dir('*_within_*');

x.a = zeros(1,length(filenames_across));
y.a = zeros(length(filenames_across),7);

x.w = zeros(1,length(filenames_within));
y.w = zeros(length(filenames_within),7);

for i = 1:length(filenames_across)
    S = load(filenames_across(i).name);
    x.a(i) = S.dataout.meandomfreq;
    y.a(i,:) = S.dataout.cidomfreqdiff_oct;
    clear S;
end
[x.a,ind] = sort(x.a);
y.a = y.a(ind,:);



for i = 1:length(filenames_within)
    S = load(filenames_within(i).name);
    x.w(i) = S.dataout.meandomfreq;
    y.w(i,:) = S.dataout.cidomfreqdiff_oct;
    clear S;
end
[x.w,ind] = sort(x.w);
y.w = y.w(ind,:);

figure;
plot(x.a,y.a,'k+','markersize',10);hold on; plot(x.w,y.w,'r+','markersize',10);hold on;

cd 'C:\LeuvenDataAnalysis';
[y.fit.a,y.std.a] = mylowessbootstrap(x.a',y.a',0.2,'lowess',200);
[y.fit.w,y.std.w] = mylowessbootstrap(x.w',y.w',0.2,'lowess',200);
