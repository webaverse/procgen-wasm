#ifndef _BIOMES_H_
#define _BIOMES_H_

#include <string>
#include <tuple>

class UV {
public:
  float u;
  float v;
};

// extern std::tuple<float, float> groundColors[];
// extern std::tuple<float, float> groundNormals[];
// extern std::tuple<float, float> groundHeights[];
// extern std::tuple<float, float> groundEmissives[];
// extern UV BIOME_UVS[];

// constexpr int waterLevel = 5;
constexpr float waterBaseHeight = 64.0f;

enum class BIOME : unsigned char {
  biOcean = 0,
  biPlains = 1,
  biDesert = 2,
  biExtremeHills = 3,
  biForest = 4,
  biTaiga = 5,
  biSwampland = 6,
  biRiver = 7,
  biNether = 8,
  biEnd = 9,
  biFrozenOcean = 10,
  biFrozenRiver = 11,
  biTundra = 12,
  biIceMountains = 13,
  biMushroomIsland = 14,
  biMushroomShore = 15,
  biBeach = 16,
  biDesertHills = 17,
  biForestHills = 18,
  biTaigaHills = 19,
  biExtremeHillsEdge = 20,
  biJungle = 21,
  biJungleHills = 22,
  biJungleEdge = 23,
  biDeepOcean = 24,
  biStoneBeach = 25,
  biColdBeach = 26,
  biBirchForest = 27,
  biBirchForestHills = 28,
  biRoofedForest = 29,
  biColdTaiga = 30,
  biColdTaigaHills = 31,
  biMegaTaiga = 32,
  biMegaTaigaHills = 33,
  biExtremeHillsPlus = 34,
  biSavanna = 35,
  biSavannaPlateau = 36,
  biMesa = 37,
  biMesaPlateauF = 38,
  biMesaPlateau = 39,
  biSunflowerPlains = 40,
  biDesertM = 41,
  biExtremeHillsM = 42,
  biFlowerForest = 43,
  biTaigaM = 44,
  biSwamplandM = 45,
  biIcePlainsSpikes = 46,
  biJungleM = 47,
  biJungleEdgeM = 48,
  biBirchForestM = 49,
  biBirchForestHillsM = 50,
  biRoofedForestM = 51,
  biColdTaigaM = 52,
  biMegaSpruceTaiga = 53,
  biMegaSpruceTaigaHills = 54,
  biExtremeHillsPlusM = 55,
  biSavannaM = 56,
  biSavannaPlateauM = 57,
  biMesaBryce = 58,
  biMesaPlateauFM = 59,
  biMesaPlateauM = 60,

  teDirt = 61,
  teStone = 62,
  
  waterRiver = 63,
  waterOcean = 64,
  waterRiverFrozen = 65,
  waterOceanFrozen = 66,

  liLava = 67,

  NUM_BIOMES = 68,
};

class Biome {
  public:
    float baseHeight;
    float amps[3][2];
    unsigned int color;
    unsigned char index;
};

const Biome BIOMES[] {
  { // biOcean
    50,
    {	   
      {
        0.1,
        2
      },
      {
        0.05,
        10
      },
      {
        0.01,
        8
      }
    },
    6316128,
    0
  },
  { // biPlains
    68,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    9286496,
    1
  },
  { // biDesert
    68,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    16421912,
    2
  },
  { // biExtremeHills
    100,
    {	   
      {
        0.2,
        4
      },
      {
        0.05,
        20
      },
      {
        0.01,
        16
      }
    },
    6316128,
    3
  },
  { // biForest
    70,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        2
      },
      {
        0.01,
        4
      }
    },
    353825,
    4
  },
  { // biTaiga
    70,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        2
      },
      {
        0.01,
        4
      }
    },
    747097,
    5
  },
  { // biSwampland
    61.5,
    {	   
      {
        0.1,
        1.1
      },
      {
        0.05,
        1.5
      },
      {
        0.02,
        2.5
      }
    },
    3145690,
    6
  },
  { // biRiver
    56,
    {	   
      {
        0.2,
        0.1
      },
      {
        0.05,
        0.1
      },
      {
        0.01,
        0.1
      }
    },
    16440917,
    7
  },
  { // biNether
    0,
    {	   
      {
        0.1,
        0
      },
      {
        0.01,
        0
      },
      {
        0.01,
        0
      }
    },
    8323072,
    8
  },
  { // biEnd
    0,
    {	   
      {
        0.1,
        0
      },
      {
        0.01,
        0
      },
      {
        0.01,
        0
      }
    },
    32767,
    9
  },
  { // biFrozenOcean
    40,
    {	   
      {
        0.1,
        2
      },
      {
        0.05,
        12
      },
      {
        0.01,
        10
      }
    },
    9474192,
    10
  },
  { // biFrozenRiver
    56,
    {	   
      {
        0.2,
        0.1
      },
      {
        0.05,
        0.1
      },
      {
        0.01,
        0.1
      }
    },
    16445632,
    11
  },
  { // biTundra
    68,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    16777215,
    12
  },
  { // biIceMountains
    80,
    {	   
      {
        0.2,
        2
      },
      {
        0.05,
        10
      },
      {
        0.01,
        8
      }
    },
    10526880,
    13
  },
  { // biMushroomIsland
    80,
    {	   
      {
        0.1,
        2
      },
      {
        0.05,
        8
      },
      {
        0.01,
        6
      }
    },
    16711935,
    14
  },
  { // biMushroomShore
    64,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        2
      },
      {
        0.01,
        4
      }
    },
    10486015,
    15
  },
  { // biBeach
    64,
    {	   
      {
        0.1,
        0.5
      },
      {
        0.05,
        1
      },
      {
        0.01,
        1
      }
    },
    16440917,
    16
  },
  { // biDesertHills
    75,
    {	   
      {
        0.2,
        2
      },
      {
        0.05,
        5
      },
      {
        0.01,
        4
      }
    },
    13786898,
    17
  },
  { // biForestHills
    80,
    {	   
      {
        0.2,
        2
      },
      {
        0.05,
        12
      },
      {
        0.01,
        10
      }
    },
    2250012,
    18
  },
  { // biTaigaHills
    80,
    {	   
      {
        0.2,
        2
      },
      {
        0.05,
        12
      },
      {
        0.01,
        10
      }
    },
    1456435,
    19
  },
  { // biExtremeHillsEdge
    80,
    {	   
      {
        0.2,
        3
      },
      {
        0.05,
        16
      },
      {
        0.01,
        12
      }
    },
    8359807,
    20
  },
  { // biJungle
    70,
    {	   
      {
        0.1,
        3
      },
      {
        0.05,
        6
      },
      {
        0.01,
        6
      }
    },
    5470985,
    21
  },
  { // biJungleHills
    80,
    {	   
      {
        0.2,
        3
      },
      {
        0.05,
        12
      },
      {
        0.01,
        10
      }
    },
    2900485,
    22
  },
  { // biJungleEdge
    70,
    {	   
      {
        0.1,
        3
      },
      {
        0.05,
        6
      },
      {
        0.01,
        6
      }
    },
    6458135,
    23
  },
  { // biDeepOcean
    40,
    {	   
      {
        0.1,
        2
      },
      {
        0.05,
        12
      },
      {
        0.01,
        10
      }
    },
    3158064,
    24
  },
  { // biStoneBeach
    40,
    {	   
      {
        0.1,
        2
      },
      {
        0.05,
        12
      },
      {
        0.01,
        10
      }
    },
    10658436,
    25
  },
  { // biColdBeach
    64,
    {	   
      {
        0.1,
        0.5
      },
      {
        0.05,
        1
      },
      {
        0.01,
        1
      }
    },
    16445632,
    26
  },
  { // biBirchForest
    70,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        2
      },
      {
        0.01,
        4
      }
    },
    3175492,
    27
  },
  { // biBirchForestHills
    80,
    {	   
      {
        0.2,
        2
      },
      {
        0.05,
        10
      },
      {
        0.01,
        8
      }
    },
    2055986,
    28
  },
  { // biRoofedForest
    70,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        2
      },
      {
        0.01,
        4
      }
    },
    4215066,
    29
  },
  { // biColdTaiga
    70,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        2
      },
      {
        0.01,
        4
      }
    },
    3233098,
    30
  },
  { // biColdTaigaHills
    80,
    {	   
      {
        0.2,
        2
      },
      {
        0.05,
        10
      },
      {
        0.01,
        8
      }
    },
    5864818,
    31
  },
  { // biMegaTaiga
    70,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        2
      },
      {
        0.01,
        4
      }
    },
    5858897,
    32
  },
  { // biMegaTaigaHills
    80,
    {	   
      {
        0.2,
        2
      },
      {
        0.05,
        10
      },
      {
        0.01,
        8
      }
    },
    5858905,
    33
  },
  { // biExtremeHillsPlus
    120,
    {	   
      {
        0.2,
        4
      },
      {
        0.05,
        20
      },
      {
        0.01,
        16
      }
    },
    5271632,
    34
  },
  { // biSavanna
    68,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    12431967,
    35
  },
  { // biSavannaPlateau
    80,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    10984804,
    36
  },
  { // biMesa
    70,
    {	   
      {
        0.2,
        2
      },
      {
        0.05,
        10
      },
      {
        0.01,
        8
      }
    },
    14238997,
    37
  },
  { // biMesaPlateauF
    80,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    11573093,
    38
  },
  { // biMesaPlateau
    80,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    13274213,
    39
  },
  { // biSunflowerPlains
    40,
    {	   
      {
        0.1,
        2
      },
      {
        0.05,
        12
      },
      {
        0.01,
        10
      }
    },
    11918216,
    40
  },
  { // biDesertM
    40,
    {	   
      {
        0.1,
        2
      },
      {
        0.05,
        12
      },
      {
        0.01,
        10
      }
    },
    16759872,
    41
  },
  { // biExtremeHillsM
    40,
    {	   
      {
        0.1,
        2
      },
      {
        0.05,
        12
      },
      {
        0.01,
        10
      }
    },
    8947848,
    42
  },
  { // biFlowerForest
    40,
    {	   
      {
        0.1,
        2
      },
      {
        0.05,
        12
      },
      {
        0.01,
        10
      }
    },
    2985545,
    43
  },
  { // biTaigaM
    40,
    {	   
      {
        0.1,
        2
      },
      {
        0.05,
        12
      },
      {
        0.01,
        10
      }
    },
    3378817,
    44
  },
  { // biSwamplandM
    60,
    {	   
      {
        1,
        3
      },
      {
        1.1,
        7
      },
      {
        0.01,
        0.01
      }
    },
    522674,
    45
  },
  { // biIcePlainsSpikes
    40,
    {	   
      {
        0.1,
        2
      },
      {
        0.05,
        12
      },
      {
        0.01,
        10
      }
    },
    11853020,
    46
  },
  { // biJungleM
    70,
    {	   
      {
        0.1,
        3
      },
      {
        0.05,
        6
      },
      {
        0.01,
        6
      }
    },
    8102705,
    47
  },
  { // biJungleEdgeM
    70,
    {	   
      {
        0.1,
        3
      },
      {
        0.05,
        6
      },
      {
        0.01,
        6
      }
    },
    6458135,
    48
  },
  { // biBirchForestM
    70,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        2
      },
      {
        0.01,
        4
      }
    },
    5807212,
    49
  },
  { // biBirchForestHillsM
    80,
    {	   
      {
        0.2,
        2
      },
      {
        0.05,
        10
      },
      {
        0.01,
        8
      }
    },
    4687706,
    50
  },
  { // biRoofedForestM
    70,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        2
      },
      {
        0.01,
        4
      }
    },
    6846786,
    51
  },
  { // biColdTaigaM
    70,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        2
      },
      {
        0.01,
        4
      }
    },
    2375478,
    52
  },
  { // biMegaSpruceTaiga
    70,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        2
      },
      {
        0.01,
        4
      }
    },
    4542270,
    53
  },
  { // biMegaSpruceTaigaHills
    80,
    {	   
      {
        0.2,
        2
      },
      {
        0.05,
        10
      },
      {
        0.01,
        8
      }
    },
    4542286,
    54
  },
  { // biExtremeHillsPlusM
    120,
    {	   
      {
        0.2,
        4
      },
      {
        0.05,
        20
      },
      {
        0.01,
        16
      }
    },
    7903352,
    55
  },
  { // biSavannaM
    68,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    15063687,
    56
  },
  { // biSavannaPlateauM
    80,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    10984820,
    57
  },
  { // biMesaBryce
    80,
    {	   
      {
        0.2,
        2
      },
      {
        0.1,
        30
      },
      {
        0.01,
        8
      }
    },
    16739645,
    58
  },
  { // biMesaPlateauFM
    80,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    14204813,
    59
  },
  { // biMesaPlateauM
    80,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    15905933,
    60
  },

  { // teDirt
    0,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    0x8d6e63,
    61
  },
  { // teStone
    0,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    0x78909c,
    62
  },
  
  { // waterRiver
    0,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    0x42a5f5,
    63
  },
  { // waterOcean
    0,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    0x42a5f5,
    64
  },
  { // waterRiverFrozen
    0,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    0x42a5f5,
    65
  },
  { // waterOceanFrozen
    0,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    0x42a5f5,
    66
  },

  { // liLava
    0,
    {	   
      {
        0.1,
        1
      },
      {
        0.05,
        1.5
      },
      {
        0.01,
        4
      }
    },
    0xef5350,
    67
  },
};

const BIOME BIOMES_TEMPERATURE_HUMIDITY[] {
		//       0                1                2                      3                      4                      5                      6                7                8                9                10                     11                    12                     13                     14                    15
		/*  0 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTundra,       BIOME::biTundra,       BIOME::biPlains,       BIOME::biPlains,       BIOME::biPlains, BIOME::biPlains, BIOME::biDesert, BIOME::biDesert, BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesert,       BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesert,
		/*  1 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTundra,       BIOME::biTundra,       BIOME::biPlains,       BIOME::biPlains,       BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biDesert, BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesert,       BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesert,
		/*  2 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTundra,       BIOME::biTundra,       BIOME::biPlains,       BIOME::biExtremeHills, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biDesert, BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesertHills,  BIOME::biDesertHills,  BIOME::biDesert,      BIOME::biDesert,
		/*  3 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTundra,       BIOME::biTundra,       BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biDesert,       BIOME::biDesert,      BIOME::biDesertHills,  BIOME::biDesertHills,  BIOME::biDesert,      BIOME::biDesert,
		/*  4 */ BIOME::biTundra, BIOME::biTundra, BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biForestHills,  BIOME::biDesert, BIOME::biDesert, BIOME::biDesert, BIOME::biDesertHills, BIOME::biDesert,
		/*  5 */ BIOME::biTundra, BIOME::biTundra, BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biPlains, BIOME::biForestHills,  BIOME::biForestHills, BIOME::biDesert, BIOME::biDesert, BIOME::biDesertHills, BIOME::biDesert,
		/*  6 */ BIOME::biTundra, BIOME::biTundra, BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biForestHills,  BIOME::biForestHills,  BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest,       BIOME::biForestHills, BIOME::biExtremeHills, BIOME::biDesert, BIOME::biDesertHills,      BIOME::biDesert,
		/*  7 */ BIOME::biTundra, BIOME::biTundra, BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biForestHills,  BIOME::biForestHills,  BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest,       BIOME::biForestHills, BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biDesertHills,      BIOME::biDesert,
		/*  8 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTaiga,        BIOME::biTaiga,        BIOME::biForest,       BIOME::biForest,       BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest,       BIOME::biForestHills, BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biPlains,      BIOME::biDesert,
		/*  9 */ BIOME::biTundra, BIOME::biTundra, BIOME::biTaiga,        BIOME::biTaiga,        BIOME::biForest,       BIOME::biForest,       BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest,       BIOME::biForestHills, BIOME::biExtremeHills, BIOME::biExtremeHills, BIOME::biPlains,      BIOME::biPlains,
		/* 10 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biTaiga,        BIOME::biIceMountains, BIOME::biForestHills,  BIOME::biForestHills,  BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
		/* 11 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biExtremeHills, BIOME::biForestHills,  BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biForest, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
		/* 12 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biIceMountains, BIOME::biIceMountains, BIOME::biExtremeHills, BIOME::biJungleHills,  BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
		/* 13 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biTaiga,        BIOME::biIceMountains, BIOME::biJungleHills,  BIOME::biJungleHills,  BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
		/* 14 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biTaiga,        BIOME::biTaiga,        BIOME::biJungle,       BIOME::biJungle,       BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
		/* 15 */ BIOME::biTaiga,  BIOME::biTaiga,  BIOME::biTaiga,        BIOME::biTaiga,        BIOME::biJungle,       BIOME::biJungle,       BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle, BIOME::biJungle,       BIOME::biJungle,      BIOME::biSwampland,    BIOME::biSwampland,    BIOME::biSwampland,   BIOME::biSwampland,
	};

// computed in the frontend
const UV BIOME_UVS[] {
  {0, 0}, // biOcean
  {0.14285714285714285, 0}, // biPlains
  {0.2857142857142857, 0}, // biDesert
  {0.42857142857142855, 0}, // biExtremeHills
  {0.5714285714285714, 0}, // biForest
  {0.7142857142857143, 0}, // biTaiga
  {0.8571428571428571, 0}, // biSwampland
  {0, 0.14285714285714285}, // biRiver
  {0.14285714285714285, 0.14285714285714285}, // biNether
  {0.2857142857142857, 0.14285714285714285}, // biEnd
  {0.42857142857142855, 0.14285714285714285}, // biFrozenOcean
  {0.5714285714285714, 0.14285714285714285}, // biFrozenRiver
  {0.7142857142857143, 0.14285714285714285}, // biTundra
  {0.8571428571428571, 0.14285714285714285}, // biIceMountains
  {0, 0.2857142857142857}, // biMushroomIsland
  {0.14285714285714285, 0.2857142857142857}, // biMushroomShore
  {0.2857142857142857, 0}, // biBeach
  {0.2857142857142857, 0.2857142857142857}, // biDesertHills
  {0.42857142857142855, 0.2857142857142857}, // biForestHills
  {0.5714285714285714, 0.2857142857142857}, // biTaigaHills
  {0.42857142857142855, 0}, // biExtremeHillsEdge
  {0.7142857142857143, 0.2857142857142857}, // biJungle
  {0.8571428571428571, 0.2857142857142857}, // biJungleHills
  {0, 0.42857142857142855}, // biJungleEdge
  {0, 0}, // biDeepOcean
  {0.14285714285714285, 0.42857142857142855}, // biStoneBeach
  {0.2857142857142857, 0.42857142857142855}, // biColdBeach
  {0.5714285714285714, 0.2857142857142857}, // biBirchForest
  {0.42857142857142855, 0.42857142857142855}, // biBirchForestHills
  {0.7142857142857143, 0}, // biRoofedForest
  {0.5714285714285714, 0.42857142857142855}, // biColdTaiga
  {0.7142857142857143, 0}, // biColdTaigaHills
  {0.7142857142857143, 0.42857142857142855}, // biMegaTaiga
  {0.8571428571428571, 0.42857142857142855}, // biMegaTaigaHills
  {0.42857142857142855, 0}, // biExtremeHillsPlus
  {0.42857142857142855, 0.42857142857142855}, // biSavanna
  {0, 0.5714285714285714}, // biSavannaPlateau
  {0.2857142857142857, 0.2857142857142857}, // biMesa
  {0.14285714285714285, 0.5714285714285714}, // biMesaPlateauF
  {0.14285714285714285, 0.5714285714285714}, // biMesaPlateau
  {0.2857142857142857, 0.5714285714285714}, // biSunflowerPlains
  {0.42857142857142855, 0.5714285714285714}, // biDesertM
  {0.42857142857142855, 0}, // biExtremeHillsM
  {0.14285714285714285, 0}, // biFlowerForest
  {0.7142857142857143, 0}, // biTaigaM
  {0.8571428571428571, 0}, // biSwamplandM
  {0.5714285714285714, 0.5714285714285714}, // biIcePlainsSpikes
  {0.7142857142857143, 0.2857142857142857}, // biJungleM
  {0, 0.42857142857142855}, // biJungleEdgeM
  {0.5714285714285714, 0.2857142857142857}, // biBirchForestM
  {0.42857142857142855, 0.42857142857142855}, // biBirchForestHillsM
  {0.7142857142857143, 0}, // biRoofedForestM
  {0.5714285714285714, 0.42857142857142855}, // biColdTaigaM
  {0.7142857142857143, 0.42857142857142855}, // biMegaSpruceTaiga
  {0.8571428571428571, 0.42857142857142855}, // biMegaSpruceTaigaHills
  {0.42857142857142855, 0}, // biExtremeHillsPlusM
  {0.42857142857142855, 0.42857142857142855}, // biSavannaM
  {0, 0.5714285714285714}, // biSavannaPlateauM
  {0.7142857142857143, 0.5714285714285714}, // biMesaBryce
  {0.14285714285714285, 0.5714285714285714}, // biMesaPlateauFM
  {0.14285714285714285, 0.5714285714285714}, // biMesaPlateauM
  {0.8571428571428571, 0.5714285714285714}, // teDirt
  {0, 0}, // teStone
  {0, 0.7142857142857143}, // liWater
  {0.14285714285714285, 0.7142857142857143}, // liWaterOcean
  {0.2857142857142857, 0.7142857142857143}, // liWaterRiverFrozen
  {0.42857142857142855, 0.7142857142857143}, // liWaterOceanFrozen
  {0.5714285714285714, 0.7142857142857143}, // liLava
};

enum class MATERIAL : uint8_t {
  GRASS = 0,
  DIRT = 1,
  STONE = 2,
  ROCK = 3,

  NUM_MATERIALS = 4,
};

bool isWaterBiome(unsigned char b);

#endif // _BIOMES_H_