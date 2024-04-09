# list of commas and the temperament family name they correspond to
# should be restricted to the 19-limit!

# fmt: off
names = [
# temps organised only by rank:
    # rank 2
    ("512/507", "512/507"),
    ("diaschismic", "2147483648/2109289329"),
    ("flattone", "137781/131072"),
    ("flattone", "729/704"),
    ("gariboh", "3125/3087"),
    ("garischismic", "[25 -14 0 -1]"),
    ("hedgehog", "118098/117649"),
    ("injera", "6561/6272"),
    ("liese", "[-9 11 0 -3]"),
    ("meantone", "59049/57344"),
    ("mohajira", "214358881/204800000"),
    ("nexus", "1771561/1769472"),
    ("octacot", "5764801/5668704"),
    ("orwell", "839808/823543"),
    ("pajara", "2197265625/1977326743"),
    ("porcupine", "1331/1296"),
    ("rastmic", "243/242"),
    ("sensamagic", "245/243"),
    ("sensi", "1647086/1594323"),
    ("superpyth", "20480/19683"),
    ("supra", "8192/8019"),
    ("tridecimal", "1053/1024"),

    # rank 3
    ("144/143", "144/143"),  # grossmic? needs better name!
    ("55/54", "55/54"), # unnamed
    ("729/728", "729/728"),  # squbemic, needs better name
    ("apollo", "100/99"), # ptolemismic
    ("arcturus", "15625/15309"),
    ("breed", "2401/2400"), # breedsmic?
    ("buzurg", "169/168"),  # aka buzurg(ism)ic, dhanvantari
    ("counterpyth", "1216/1215"),
    ("demeter", "686/675"), # sengic
    ("guanyin", "1728/1715"),  # orwellismic
    ("indra", "16875/16807"), # mirkwai
    ("island", "676/675"),
    ("konbini", "56/55"),
    ("landscape", "250047/250000"),
    ("parapyth", "352/351"),
    ("marvel", "225/224"),
    ("marveltwin", "325/324"),
    ("minerva", "5632/5625"), # vishdel
    ("mothwellsmic", "99/98"),  # needs better name!
    ("odin", "[-17 24 -18 0 6]"),
    ("pele", "5120/5103"),  # hemifamity
    ("akea", "2200/2187"),  # small tetracot diesis
    ("pentacircle", "896/891"),
    ("pine", "4000/3993"),
    ("ragismic", "4375/4374"),
    ("starling", "126/125"),
    ("supermagic", "875/864"),
    ("symbiotic", "19712/19683"),
    ("thor", "1890625/1889568"),
    ("trimitone", "8019/8000"),
    ("vulkan", "512/495"),
    ("zeus", "121/120"),  # biyatismic
    ("zeus", "176/175"),  # valinorismic
    ("zeus", "6144/6125"),  # porwell

    # rank 4
    ("animist", "105/104"),
    ("biome", "91/90"),
    ("gentle", "364/363"), # minor minthmic
    ("huntmic", "640/637"), # ehh
    ("kalismic", "9801/9800"),
    ("keenanismic", "385/384"), # keenan? (i can literally ask him) "undecimal kleismic" doesnt sound that nice
    ("lehmerismic", "3025/3024"), # lehmer?
    ("mynucumic", "196/195"), # needs better name!
    ("ratwolf", "351/350"),
    ("swetismic", "540/539"), # swets?
    ("werckismic", "441/440"), # needs better name!


# organised temps below:

# continuums implicitly also included (that i"m aware of):
#
# schismic-pythagorean: schisma^n = pyth comma (high accuracy):
# gracecordial (-1), compton (0), meantone (1), diaschismic (2), misty (3), undim (4) [...] schismic (infinity)
# chromatic-diatonic: (25/24)^n = 16/15 (v low accuracy except wurschmidt):
# father (0), augmented (1), wurschmidt (1.5), magic (2) [...] dicot (infinity)

# how many 2187/2048"s (pyth chromatic semitones) are in a 256/243 (pyth diatonic semitone)?
# 2.3 subgroup: (this continuum is significant as it encompasses all 5L2s tunings) 
("blackwood","256/243"),# 0 (5 EDO)
("gothic","[27 -17]"),# 1/2 (17 EDO)
("mystery","[46 -29]"),# 2/3 (29 EDO)
("countercomp","[65 -41]"),# 3/4 (41 EDO)
("mercator","[-84 53]"), # 4/5 (53 EDO)
("compton","531441/524288"),# 1 (12 EDO)
("19-comma","[-30 19]"),# 2 AKA 19 & 19c; in the syntonic-kleismic continuum at n=0 & supported up to 76 EDO by pval
("whitewood","2187/2048"),# infinity (7 EDO)

# 2.3.5 subgroup:
("father","16/15"),# for chromatic-diatonic continuum"s completeness
# syntonic-chromatic: how many 81/80"s are in a 25/24? (rank 2 except whitewood (see above) at -2)
("mavila","135/128"),# -1
("dicot","25/24"),# 0
("porcupine","250/243"),# 1
("tetracot","20000/19683"),# 2
("artoneutral","32000000000/31381059609"),# 2.5; significant as it"s the 87&94 temp
("amity","1600000/1594323"),# 3
# ideal number of 81/80"s is about here
("undetrita","205891132094649/204800000000000"),# 3.5; significant as it"s the 111&118 temp
("gravity","129140163/128000000"),# 4
("absurdity", "10460353203/10240000000"),# 5
# how many 81/80"s in 128/125?
("meantone","81/80"),# infinity
("python","43046721/41943040"), # 4
# compton is 3
("gracecordial","[-34 20 1]"),# 2.5
("schismic","32805/32768"),# 2
# ideal number of 81/80"s is ~1.9
("undim","[41 -20 -4]"),# 5/3
("misty","67108864/66430125"),# 1.5
("diaschismic","2048/2025"), # 1
# augmented is 0; see 2.5 subgroup
("diminished","648/625"),# -1
# syntonic-kleismic: how many 81/80"s in the 19-comma? (the 19-comma is 0)
("negri","16875/16384"),# 4
("magic","3125/3072"),# 5
("kleismic","15625/15552"),# 6
("enneadecal","[-14 -19 19]"),# 6+1/3; 19 = 28/27
# ideal number of 81/80"s is about here
("parakleismic","[8 14 -13]"),# 6+1/2
("sensi","78732/78125"),# 7
("unicorn","1594323/1562500"),# 8
# syntonic-31: how many 81/80"s in the 31-comma ([-49 31])?
# meantone is infinity
("nusecond","[5 13 -11]"),# 11
("myna","[9 9 -10]"),# 10
("valentine","[13 5 -9]"),# 9
("würschmidt","[17 1 -8]"),# 8
# birds is here at 7.75 btw but it"s in the no-3"s 5-limit
("hemithirds","[38 -2 -15]"),# 7.5 (AKA luna; very close to just)
("tertiaseptal","[-59 5 22]"),# 7.333... = 22/3
("orwell","[-21 3 7]"),# 7 (AKA orson; corresponds to the (strangely-named) "semicomma")
("miracle","[-25 7 6]"),# 6 (AKA ampersand)
("tritonic","[-29 11 5]"),# 5
("sentinel","[-33 15 4]"),# 4 (as we want to include nusecond(?) may as well include this too as it"s supported by more EDOs)
# period- or gen-based microtemps (enneadecal is above):
("ennealimmal","[1 -27 18]"),# 1\9 = 27/25; note: enneadeca / ennealimma = schisma
("chlorine","[-52 -17 34]"),# 1 = 25/24
("vishnu","[23 6 -14]"),# (4/3)/(25/24)^7

# 2.5 subgroup:
("bug","27/25"),
("augmented","128/125"), # 1
("birds","[72 0 -31]"), # 10

# 2.7 subgroup:
("cloudy","[-14 0 0 5]"), # 4
("birds","[-87 0 0 31]"), # 25

# 2.5.7 subgroup:
# if u want u can think of rainy as part of the continuum from cloudy (one extreme) to augmented (the other)
# but it"s basically trivial as 128/125 = [-14 0 0 5] is very near 1 so other tunings tend to be a bit silly
("rainy","2100875/2097152"),#   128/125 = [-14 0 0 5] (cloudly comma)
("miracle","823543/819200"),# 2.5.7[10&31]-comma
("didacus","3136/3125"),# = (50/49)/(128/125)
("jubilismic","50/49"),

# 2.3.7 subgroup:
# archytas/septimal-diatonic: how many 64/63"s in 256/243?
# blackwood is 0; see 2.3 subgroup
("trienstonian","28/27"),# 1 (included for completeness)
("semaphore","49/48"),# 2
("slendric","1029/1024"),# 3; maybe a good argument for renaming to "gamelismic" (contrast "slendrismic")
("slendrismic","[36 -5 0 -10]"),# 3.33... = 10/3; epic temp making 147/128 = 1; near-just (~3.3)
("septiness","67108864/66706983"),# 3.5
("buzzard","65536/64827"),# 4
("superpyth","64/63"),# infinity
# period-based microtemps (think these are nowhere in the continuum above?)
("ennealimmal","[-11 -9 0 9]"), # 2\9 = 7/6
("enneadecal","[-37 57 0 -19]"), # 19 = 28/27

]
