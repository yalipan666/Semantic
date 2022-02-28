# -*- coding: utf-8 -*-
"""
Created on Wed Feb 23 16:22:36 2022

@author: panyz
"""

import seaborn as sns
sns.set(style = "ticks",palette="muted", font='arial', color_codes=True)
sns.axes_style("ticks", {"xtick.major.size": 7, "ytick.major.size": 7})
import matplotlib.pyplot as plt


x = [-0.00197108567966008,	-0.00467241687901599,	-0.0142202766309802,	-0.00288174791202399,	0.00356254445714156,	-0.00275006894972543,	-0.00291793217200571,	-0.00318348276520496,	-0.00472156170574349,	-0.0175787175436676,	0.000283286656495275,	0.00129944917323996,	0.0124892138117057,	0.000532439889682095,	0.00720350467294630,	-0.000113940275377381,	-0.000616486613489176,	-0.0106226614831061,	-0.00197636398618489,	-0.000939961081150423,	0.00386980659353443,	0.00281208371659599,	0.00417181271978365,	-0.0144263645628962,	-0.00599376605761048,	-0.0106346284091592,	-0.0110291728711136,	-0.0221091326461002,	-0.00791996210985320]
y = [-2.23938898087018,	-2.34865675876022,	-2.27178760511396,	-3.20430074200316,	-1.59235796225801,	-1.13692121795389,	-2.25868890118371,	-1.18717244768286,	-0.541641410309349,	-2.59051345360796,	-2.98502968744687,	-5.13358656382997,	-2.03932940270958,	-2.18945296406217,	-2.90079979252592,	0.972941165937062,	-2.64494291570624,	-1.29328667762229,	-0.848486633292509,	0.0705304695243765,	-3.99498833124664,	-1.62936242705460,	-0.573019828247738,	-2.07605183494035,	-2.89894869722889,	-0.646967429554431,	-1.10588954281125,	-0.648110794934207,	-4.03548876090234]
sns.jointplot(x, y, kind="reg", truncate=False, color="grey", height = 3)
plt.savefig('U:/writing/Parafoveal_Semantic/figs/TagPre_N400Targ.svg', dpi = 300)


y = [229.379043109668,	235.218033633034,	154.583788780664,	238.941246843434,	231.767757190726,	201.665509594572,	180.824965850122,	197.517988157676,	123.348585875930,	290.998034586941,	155.674896527084,	232.989393453768,	222.215714736652,	155.972927888084,	154.594683042652,	172.033654331779,	308.966580312049,	287.369936417749,	250.411944721944,	162.245621600622,	212.684895694583,	224.800800172050,	224.034691298285,	276.872573606949,	272.346896957209,	331.145371121934,	280.965885555417,	341.686756750194,	315.797696643634]
sns.jointplot(x, y, kind="reg", truncate=False, color="grey", height = 3)
plt.savefig('U:/writing/Parafoveal_Semantic/figs/TagPre_DurationPerwrd.svg', dpi = 300)


y = [37.2927180966114,	25.0794912559618,	21.3114754098361,	38.4615384615385,	33.8675213675214,	26.9754117927906,	33.7719298245614,	22.8016749143510,	9.72222222222222,	68.7565858798736,	15.6250000000000,	31.0239817282071,	23.3766233766234,	42.9787234042553,	17.1095571095571,	3.26656394453005,	63.2119682768565,	36.2500000000000,	31.9444444444444,	15.7545271629779,	66.3278750235272,	28.7671232876712,	18.1326604181687,	18.7601428107757,	47.7280103862382,	14.9452736318408,	58.2991620609734,	9.09090909090909,	55.1282051282051]
sns.jointplot(x, y, kind="reg", truncate=False, color="grey", height = 3)
plt.savefig('U:/writing/Parafoveal_Semantic/figs/TagPre_RegdifTarg.svg', dpi = 300)


x = [-2.23938898087018,	-2.34865675876022,	-0.445309307599853,	-2.27178760511396,	-3.20430074200316,	-1.59235796225801,	-1.13692121795389,	-2.25868890118371,	-1.18717244768286,	-0.541641410309349,	-2.59051345360796,	-2.98502968744687,	-5.42713737996186,	-3.62527357099112,	-5.13358656382997,	-1.75798391313867,	-2.03932940270958,	-2.18945296406217,	-2.90079979252592,	0.972941165937062,	-1.37539616138493,	-2.64494291570624,	-1.29328667762229,	-0.848486633292509,	0.0705304695243765,	-3.99498833124664,	-1.62936242705460,	-0.573019828247738,	-2.07605183494035,	-2.89894869722889,	-0.646967429554431,	-1.10588954281125,	-0.648110794934207,	-4.03548876090234]
y = [229.379043109668,	235.218033633034,	229.548291812354,	154.583788780664,	238.941246843434,	231.767757190726,	201.665509594572,	180.824965850122,	197.517988157676,	123.348585875930,	290.998034586941,	155.674896527084,	227.283695471196,	272.622995580808,	232.989393453768,	211.124693969225,	222.215714736652,	155.972927888084,	154.594683042652,	172.033654331779,	208.426880896881,	308.966580312049,	287.369936417749,	250.411944721944,	162.245621600622,	212.684895694583,	224.800800172050,	224.034691298285,	276.872573606949,	272.346896957209,	331.145371121934,	280.965885555417,	341.686756750194,	315.797696643634]
sns.jointplot(x, y, kind="reg", truncate=False, color="grey", height = 3)
plt.savefig('U:/writing/Parafoveal_Semantic/figs/N400Targ_DurationPerwrd.svg', dpi = 300)


y = [37.2927180966114,	25.0794912559618,	8.86075949367088,	21.3114754098361,	38.4615384615385,	33.8675213675214,	26.9754117927906,	33.7719298245614,	22.8016749143510,	9.72222222222222,	68.7565858798736,	15.6250000000000,	37.0129870129870,	27.3972602739726,	31.0239817282071,	7.16666666666667,	23.3766233766234,	42.9787234042553,	17.1095571095571,	3.26656394453005,	26.9041769041769,	63.2119682768565,	36.2500000000000,	31.9444444444444,	15.7545271629779,	66.3278750235272,	28.7671232876712,	18.1326604181687,	18.7601428107757,	47.7280103862382,	14.9452736318408,	58.2991620609734,	9.09090909090909,	55.1282051282051]
sns.jointplot(x, y, kind="reg", truncate=False, color="grey", height = 3)
plt.savefig('U:/writing/Parafoveal_Semantic/figs/N400Targ_RegdifTarg.svg', dpi = 300)







