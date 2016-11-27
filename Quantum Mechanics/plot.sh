cd Well
mkdir Images

gnuplot << abc
set terminal png
set xr [0:16000]
set yr [0:0.004]

list = system('ls schrodinger*')
i = 1

do for [file in list] {
    set output sprintf('Images/%s.png', file)
    set title sprintf("%s", file, i)
    plot file with lines
    i = i + 1
}

abc

cd ..
ffmpeg -i Well/Images/schrodinger%01d.png output.gif
rm -r Well/Images