from typing import TypeAlias

TempDict: TypeAlias = dict[tuple[int, ...], str]
names: dict[tuple[int, ...], TempDict] = {}

# fmt: off

# rank 1

names[(2, 3)] = {
    (5, 8): "blackwood",
    (7, 11): "whitewood",
    (12, 19): "compton",
}

names[(2, 5)] = {
    (3, 7): "augmented",
}

names[(2, 7)] = {
    (5, 14): "cloudy",
}

names[(3, 5)] = {
    (2, 3): "bug",
}

# rank 2

names[(2, 3, 5)] = {
    (1, 0, -4, 0, 1, 4):   "meantone",
    (1, 0, 1, 0, 6, 5):    "kleismic",  # hanson?
    (1, 0, 15, 0, 1, -8):  "schismatic",
    (1, 0, 2, 0, 5, 1):    "magic",
    (1, 0, 23, 0, 1, -13): "superpyth",
    (1, 0, 3, 0, 7, -3):   "orwell",
    (1, 0, 4, 0, 1, -1):   "father",
    (1, 0, 7, 0, 1, -3):   "mavila",
    (1, 1, 0, 0, 2, 8):    "mohajira",
    (1, 1, 1, 0, 4, 9):    "tetracot",
    (1, 1, 1, 0, 8, 18):   "octacot",
    (1, 1, 2, 0, 2, 1):    "dicot",
    (1, 1, 2, 0, 9, 5):    "valentine",
    (1, 1, 3, 0, 6, -7):   "miracle",
    (1, 2, 2, 0, 4, -3):   "negri",
    (1, 2, 3, 0, 3, 5):    "porcupine",
    (1, 6, 8, 0, 7, 9):    "sensi",
    (1, 7, 3, 0, 8, 1):    "würschmidt",
    (1, 9, 9, 0, 10, 9):   "myna",
    (2, 0, -8, 0, 1, 4):   "injera",
    (2, 0, 11, 0, 1, -2):  "diaschismic",  # aka pajara
    (2, 1, 1, 0, 3, 5):    "hedgehog",
    (4, 0, 3, 0, 1, 1):    "diminished",
    (9, 1, 1, 0, 2, 3):    "ennealimmal",
}

names[(2, 3, 7)] = {
    (1, 0, -13, 0, 1, 10): "meantone",
    (1, 0, 1, 0, 7, 8):    "orwell",
    (1, 0, 17, 0, 1, -9):  "flattone",
    (1, 0, 2, 0, 2, 1):    "semaphore",
    (1, 0, 6, 0, 1, -2):   "superpyth",
    (1, 1, 3, 0, 3, -1):   "slendric",
    (1, 1, 3, 0, 6, -2):   "miracle",
    (1, 6, 11, 0, 7, 13):  "sensi",
    (2, 0, 12, 0, 1, -2):  "pajara",
    (2, 0, 31, 0, 1, -8):  "diaschismic",
    (9, 1, 12, 0, 2, 2):   "ennealimmal",
}

names[(2, 3, 11)] = {
    (1, 0, -6, 0, 1, 6):  "flattone",
    (1, 0, 13, 0, 1, -6): "supra", # superpyth adjacent, 8192/8019
    (1, 1, 2, 0, 2, 5):   "rastmic",
    (1, 2, 4, 0, 3, 4):   "porcupine",
}

names[(2, 3, 13)] = {
    (1, 0, 10, 0, 1, -4): "tridecimal", # 1053/1024
    (1, 1, 4, 0, 2, -1):  "512/507-comma", # unnamed
}


names[(2, 5, 7)] = {
    (1, 0, -3, 0, 2, 5): "didacus",  # hemimean?
    (1, 2, 3, 0, 5, -3): "rainy",
    (1, 3, 3, 0, 7, 2):  "miracle",  # quince / mercy
    (2, 0, 1, 0, 1, 1):  "jubilismic",
}

names[(3, 5, 7)] = {
    (1, 1, 1, 0, 3, 5):  "gariboh",  # what?
    (1, 1, 2, 0, 2, -1): "sensamagic",  # Bohlen-Pierce? Lambda?
}

# rank 3

names[(2, 3, 5, 7)] = {
    (1, 0, 0, -1, 0, 1, 0, -2, 0, 0, 1, 3): "starling",
    (1, 0, 0, -5, 0, 1, 0, 2, 0, 0, 1, 2):  "marvel",
    (1, 0, 0, 1, 0, 1, 0, 7, 0, 0, 1, -4):  "landscape",
    (1, 0, 0, 10, 0, 1, 0, -6, 0, 0, 1, 1): "hemifamity",  # this name name is so awful
    (1, 0, 0, 2, 0, 1, 0, 1, 0, 0, 3, -1):  "orwellismic", # bad. confusing wrt orwell
    (1, 0, 0, 5, 0, 1, 0, 3, 0, 0, 1, -3):  "supermagic", # keemic?
    (1, 0, 1, 4, 0, 1, 1, -1, 0, 0, 2, -3): "porwell", # lame, zeus would be cooler
    (1, 1, 1, 2, 0, 2, 1, 1, 0, 0, 2, 1):   "breed",  # jove? breedsmic (ew)?
}

names[(2, 3, 5, 11)] = {
    (1, 0, 0, 1, 0, 1, 0, 3, 0, 0, 1, -1): "55/54-comma",  # unnamed, 'telepathmic'??
    (1, 0, 0, 2, 0, 1, 0, -2, 0, 0, 1, 2): "ptolemismic",  # why not just ptolemaic
    (1, 0, 0, 6, 0, 1, 0, -6, 0, 0, 1, 3): "trimitone",
    (1, 0, 1, 2, 0, 1, 1, 1, 0, 0, 2, 1):  "biyatismic",
    (1, 2, 0, 1, 0, 3, 0, -1, 0, 0, 1, 1): "pine",
}

names[(2, 3, 7, 11)] = {
    (1, 0, 0, 1, 0, 1, 0, -2, 0, 0, 1, 2): "mothwellsmic", # aggressive eye-roll, can we get something better for my dear 99/98
    (1, 0, 0, 7, 0, 1, 0, -4, 0, 0, 1, 1): "pentacircle",

}

names[(2, 5, 7, 11)] = {
    (1, 0, 0, -4, 0, 1, 0, 2, 0, 0, 1, 1): "valinorismic", # do better. shows up a lot!
    (1, 0, 0, 3, 0, 1, 0, -1, 0, 0, 1, 1): "konbini",
}

names[(2, 3, 5, 13)] = {
    (1, 0, 0, -1, 0, 2, 0, 3, 0, 0, 1, 1): "island",
    (1, 0, 0, 2, 0, 1, 0, 4, 0, 0, 1, -2): "marveltwin",
}

names[(2, 3, 11, 13)] = {
    (1, 0, 0, 4, 0, 1, 0, 2, 0, 0, 1, -1): "grossmic", # comma is named grossma, so thats how it is
    (1, 0, 0, 5, 0, 1, 0, -3, 0, 0, 1, 1): "major minthmic", # id prefer something simpler for this
}


names[(2, 3, 7, 13)] = {
    (1, 0, 1, 2, 0, 1, 1, 1, 0, 0, 2, 1):   "buzurgic",
    (1, 0, 0, -3, 0, 1, 0, 6, 0, 0, 1, -1): "squbema", # ???
}


# rank 4

names[(2, 3, 5, 7, 11)] = {
    (1, 0, 0, 0, -3, 0, 1, 0, 0, 2, 0, 0, 1, 0, -1, 0, 0, 0, 1, 2): "werckismic",
    (1, 0, 0, 0, 2, 0, 1, 0, 0, 3, 0, 0, 1, 0, 1, 0, 0, 0, 1, -2):  "swetismic", # lmao
    (1, 0, 0, 0, 2, 0, 1, 0, 1, 2, 0, 0, 1, 0, -1, 0, 0, 0, 2, 1):  "lehmerismic", # why not just lehmer
    (1, 0, 0, 0, 7, 0, 1, 0, 0, 1, 0, 0, 1, 0, -1, 0, 0, 0, 1, -1): "keenanismic",
    (2, 0, 0, 0, 3, 0, 1, 0, 0, -2, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1):  "kalismic",
}

names[(2, 3, 5, 7, 13)] = {
    (1, 0, 0, 0, -3, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1):  "animist",
    (1, 0, 0, 0, 1, 0, 1, 0, 0, -3, 0, 0, 1, 0, 2, 0, 0, 0, 1, 1):  "ratwolf",
    (1, 0, 0, 0, 1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0, 0, 0, 1, -1):  "biome",
    (1, 0, 0, 0, 2, 0, 1, 0, 0, -1, 0, 0, 1, 0, -1, 0, 0, 0, 1, 2): "mynucumic", # ???
    (1, 0, 0, 0, 7, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, -2):  "huntmic",
}

names[(2, 3, 7, 11, 13)] = {
    (1, 0, 0, 0, -2, 0, 1, 0, 0, 1, 0, 0, 1, 0, -1, 0, 0, 0, 1, 2): "minor minthmic",
}
