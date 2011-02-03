# File name: saveplot - creates a PostScript file using  Gnuplot version 3.8
# to save the current plot as a postscript file issue the commands:
#  gnuplot>   set out 'plotfile.ps'
#  gnuplot>   load 'saveplot'
set size 1.0 , 0.6
set terminal postscript portrait enhanced color solid lw 2 "Helvetica" 14 
set output "my-plot.ps"
replot
set terminal x11
set size 1,1

#      set terminal postscript {<mode>} {enhanced | noenhanced}
#                              {color | colour | monochrome}
#                              {blacktext | colortext | colourtext}
#                              {solid | dashed} {dashlength | dl <DL>}
#                              {linewidth | lw <LW>}
#                              {<duplexing>}
#                              {"<fontname>"} {<fontsize>}


#     set terminal gif {transparent} {interlace}
#                       {tiny | small | medium | large | giant}
#                       {size <x>,<y>}
#                       {<color0> <color1> <color2> ...}

#     set terminal png
#             {{no}transparent} {{no}interlace}
#             {tiny | small | medium | large | giant}
#             {font <face> {<pointsize>}}
#             {size <x>,<y>} {{no}crop}
#             {{no}enhanced}
#             {<color0> <color1> <color2> ...}


