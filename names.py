from typing import TypeAlias

TempDict: TypeAlias = dict[tuple[int, ...], str]
names: dict[tuple[int, ...], TempDict] = {}

# fmt: off

# rank 1

names[(2, 3)] = {
    (5, 8): "blackwood", # 256/243 (3\5)
    (7, 11): "whitewood", # 2187/2048 (4\7)
    (12, 19): "compton", # 531441/524288 (7\12)
}

names[(2, 5)] = {
    (3, 7): "augmented", # 128/125 (1\3)
}

names[(2, 7)] = {
    (5, 14): "cloudy", # 16807/16384 (4\5)
}

names[(3, 5)] = {
    (2, 3): "bug", # 27/25
}

# rank 2

names[(2, 3, 5)] = {
    (1, 0, -4, 0, 1, 4): "meantone", # 81/80
    (1, 0, 1, 0, 6, 5): "kleismic",  # 15625/15552
    (1, 0, 15, 0, 1, -8): "schismatic", # 32805/32768
    (1, 0, 2, 0, 5, 1): "magic", # 3125/3072
    (1, 0, 23, 0, 1, -13): "superpyth", # 20480/19683
    (1, 0, 3, 0, 7, -3): "orwell", # 2109375/2097152 (semicomma)
    (1, 0, 4, 0, 1, -1): "father", # 16/15
    (1, 0, 7, 0, 1, -3): "mavila", # 135/128
    (1, 1, 0, 0, 2, 8): "mohajira",
    (1, 1, 1, 0, 4, 9): "tetracot", # 20000/19683
    (1, 1, 1, 0, 8, 18): "octacot",
    (1, 1, 2, 0, 2, 1): "dicot", # 25/24
    (1, 1, 2, 0, 9, 5): "valentine", # (6/5)/(25/24)^4 = 1990656/1953125
    (1, 1, 3, 0, 6, -7): "miracle", # 34171875/33554432 = (3/2)/(16/15)^6 (ampersand)
    (1, 2, 2, 0, 4, -3): "negri", # 16875/16384
    (1, 2, 3, 0, 3, 5): "porcupine", # 250/243
    (1, 6, 8, 0, 7, 9): "sensipent", # 78732/78125
    (1, 7, 3, 0, 8, 1): "würschmidt", # 393216/390625
    (1, 9, 9, 0, 10, 9): "myna", # (6/5)^8 / (25/6) = 10077696/9765625
    (2, 0, -8, 0, 1, 4): "injera",
    (2, 0, 11, 0, 1, -2): "diaschismic",  # aka pajara
    (2, 1, 1, 0, 3, 5): "hedgehog",
    (4, 0, 3, 0, 1, 1): "diminished", # 648/625
    (9, 1, 1, 0, 2, 3): "ennealimmal", # 27/25 = 1\9
}

names[(2, 3, 7)] = {
    (1, 0, -13, 0, 1, 10): "meantone", # 81/80
    (1, 0, 1, 0, 7, 8): "orwell", # (3/1)/(7/6)^7 = 839808/823543
    (1, 0, 17, 0, 1, -9): "flattone",
    (1, 0, 2, 0, 2, 1): "semaphore", # 49/48
    (1, 0, 6, 0, 1, -2): "superpyth", # 64/63
    (1, 1, 3, 0, 3, -1): "slendric", # 1029/1024
    (1, 1, 3, 0, 6, -2): "miracle",
    (1, 6, 11, 0, 7, 13): "sensi",
    (2, 0, 12, 0, 1, -2): "pajara",
    (2, 0, 31, 0, 1, -8): "diaschismic", # 2048/2025
    (9, 1, 12, 0, 2, 2): "ennealimmal", # 2\9 = 7/6
}

names[(2, 3, 11)] = {
    (1, 1, 2, 0, 2, 5): "rastmic", # 243/242
    (1, 0, -6, 0, 1, 6): "flattone",
    (1, 2, 4, 0, 3, 4): "porcupine", # (4/3) / (12/11)^3
    (1, 0, 13, 0, 1, -6): "supra" # superpyth adjacent, 8192/8019
}

names[(2, 3, 13)] = {
    (1, 1, 4, 0, 2, -1): "512/507-comma", # unnamed
    (1, 0, 10, 0, 1, -4): "tridecimal" # 1053/1024
}


names[(2, 5, 7)] = {
    (1, 0, -3, 0, 2, 5): "didacus", # 3136/3125; hemimean?
    (1, 2, 3, 0, 5, -3): "rainy", # 2100875/2097152
    (1, 3, 3, 0, 7, 2): "miracle", # quince / mercy
    (2, 0, 1, 0, 1, 1): "jubilismic", # 50/49
}

names[(3, 5, 7)] = {
    (1, 1, 2, 0, 2, -1): "sensamagic", # 245/243; Bohlen-Pierce? Lambda?
    (1, 1, 1, 0, 3, 5): "gariboh",  # 3125/3087; what?
}

# rank 3

names[(2, 3, 5, 7)] = {
    (1, 0, 0, -1, 0, 1, 0, -2, 0, 0, 1, 3): "starling", # 126/125
    (1, 0, 0, -5, 0, 1, 0, 2, 0, 0, 1, 2): "marvel", # 225/224
    (1, 0, 0, 1, 0, 1, 0, 7, 0, 0, 1, -4): "landscape", # 250047/250000
    (1, 0, 0, 10, 0, 1, 0, -6, 0, 0, 1, 1): "hemifamity",  # 5120/5103 this name name is so awful
    (1, 0, 0, 2, 0, 1, 0, 1, 0, 0, 3, -1): "orwellismic", # 1728/1715
    (1, 0, 0, 5, 0, 1, 0, 3, 0, 0, 1, -3): "supermagic", # 875/864, keemic?
    (1, 0, 1, 4, 0, 1, 1, -1, 0, 0, 2, -3): "porwell", # lame, zeus would be cooler
    (1, 1, 2, 5, 0, 1, 0, -1, 0, 0, 1, -5): "horwell", # 65625/65536
    (1, 1, 1, 2, 0, 2, 1, 1, 0, 0, 2, 1): "breed",  # jove? breedsmic (ew)?
    (1, 1, 2, 0, 0, 1, 0, 7, 0, 0, 1, -4): "ragismic" # 4375/4374
}

names[(2, 3, 5, 11)] = {
    (1, 0, 0, 1, 0, 1, 0, 3, 0, 0, 1, -1): "55/54-comma", # afaik no name established, 'telepathmic'??
    (1, 0, 0, 2, 0, 1, 0, -2, 0, 0, 1, 2): "ptolemismic", # 100/99; why not just ptolemaic
    (1, 0, 1, 2, 0, 1, 1, 1, 0, 0, 2, 1): "biyatismic", # 121/120
    (1, 2, 0, 1, 0, 3, 0, -1, 0, 0, 1, 1): "pine", # 4000/3993
}

names[(2, 3, 7, 11)] = {
    (1, 0, 0, 1, 0, 1, 0, -2, 0, 0, 1, 2): "mothwellsmic", # aggressive eye-roll, can we get something better for my dear 99/98
    (1, 0, 0, 7, 0, 1, 0, -4, 0, 0, 1, 1): "pentacircle", # 896/891

}

names[(2, 5, 7, 11)] = {
    (1, 0, 0, -4, 0, 1, 0, 2, 0, 0, 1, 1): "valinorismic", # do better. shows up a lot!
    (1, 0, 0, 3, 0, 1, 0, -1, 0, 0, 1, 1): "konbini",
}

names[(2, 3, 5, 13)] = {
    (1, 0, 0, 2, 0, 1, 0, 4, 0, 0, 1, -2): "marveltwin", # 325/324
    (1, 0, 0, -1, 0, 2, 0, 3, 0, 0, 1, 1): "island", # 676/675
}

names[(2, 3, 11, 13)] = {
    (1, 0, 0, 5, 0, 1, 0, -3, 0, 0, 1, 1): "major minthmic", # id prefer something simpler for this
}


# TODO: rank 4
