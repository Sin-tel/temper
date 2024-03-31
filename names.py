from typing import TypeAlias

names: dict[tuple[int, ...], dict[tuple[int, ...], str]] = {}

names[(2, 3)] = {
    (5, 8): "blackwood",
    (7, 11): "whitewood",
    (12, 19): "compton",
}

names[(2, 3, 5)] = {
    (1, 0, -4, 0, 1, 4): "meantone",
    (1, 0, 0, 0, 2, 3): "bug",
    (1, 0, 1, 0, 6, 5): "kleismic",  # hanson?
    (1, 0, 15, 0, 1, -8): "schismatic",
    (1, 0, 2, 0, 5, 1): "magic",
    (1, 0, 23, 0, 1, -13): "superpyth",
    (1, 0, 3, 0, 7, -3): "orwell",
    (1, 0, 4, 0, 1, -1): "father",
    (1, 0, 7, 0, 1, -3): "mavila",
    (1, 1, 0, 0, 2, 8): "mohajira",
    (1, 1, 1, 0, 4, 9): "tetracot",
    (1, 1, 1, 0, 8, 18): "octacot",
    (1, 1, 2, 0, 2, 1): "dicot",
    (1, 1, 3, 0, 6, -7): "miracle",
    (1, 2, 3, 0, 3, 5): "porcupine",
    (1, 6, 8, 0, 7, 9): "sensi",
    (1, 7, 3, 0, 8, 1): "würschmidt",
    (2, 0, 11, 0, 1, -2): "pajara",  # aka diaschismatic, diaschismic, srutal
    (2, 1, 1, 0, 3, 5): "hedgehog",
    (3, 0, 7, 0, 1, 0): "augmented",
    (4, 0, 3, 0, 1, 1): "diminished",
    (9, 1, 1, 0, 2, 3): "ennealimmal",
}

names[(2, 3, 7)] = {
    (1, 0, -13, 0, 1, 10): "meantone",
    (1, 0, 1, 0, 7, 8): "orwell",
    (1, 0, 17, 0, 1, -9): "flattone",
    (1, 0, 2, 0, 2, 1): "semaphore",
    (1, 0, 6, 0, 1, -2): "superpyth",
    (1, 1, 3, 0, 3, -1): "slendric",
    (1, 1, 3, 0, 6, -2): "miracle",
    (1, 6, 11, 0, 7, 13): "sensi",
    (9, 1, 12, 0, 2, 2): "ennealimmal",
}

names[(2, 3, 11)] = {
    (1, 1, 2, 0, 2, 5): "rastmic",
    (1, 0, -6, 0, 1, 6): "flattone",
}

names[(2, 5, 7)] = {
    (1, 0, -3, 0, 2, 5): "jubilismic",
}

names[(3, 5, 7)] = {
    (1, 1, 2, 0, 2, -1): "sensamagic",  # Bohlen-Pierce? Lambda?
}

names[(2, 3, 5, 7)] = {
    (1, 0, 0, -5, 0, 1, 0, 2, 0, 0, 1, 2): "marvel",
    (1, 0, 0, -1, 0, 1, 0, -2, 0, 0, 1, 3): "starling",
}

names[(2, 3, 11, 13)] = {
    (1, 0, 0, 5, 0, 1, 0, -3, 0, 0, 1, 1): "major minthmic",
}
