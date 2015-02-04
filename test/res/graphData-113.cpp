// example graph data for 113 procs
struct exGraphData {
static int exVertglbtabRaw[114];

static int exVertloctabRaw[245];
static int exVertwgtcmpRaw[245];
static int exVertwgtvisRaw[245];

static int exEdgeloctabRaw[482];
static int exEdgewgtcmpRaw[482];
static int exEdgewgtvisRaw[482];
};

int exGraphData::exVertglbtabRaw[114] = {   0,    1,    2,    3,    4,    5,    7,    8,    9,   10,   11,   12,   14,   15,   16,   17,   18,   20,   21,   22,   23,   24,   25,   27,   28,   29,   30,   31,   32,   34,   35,   36,   37,   38,   40,   41,   42,   43,   44,   45,   47,   48,   49,   50,   51,   52,   54,   55,   56,   57,   58,   60,   61,   62,   63,   64,   65,   67,   68,   69,   70,   71,   72,   74,   75,   76,   77,   78,   80,   81,   82,   83,   84,   85,   87,   88,   89,   90,   91,   92,   94,   95,   96,   97,   98,  100,  101,  102,  103,  104,  105,  107,  108,  109,  110,  111,  112,  114,  115,  116,  117,  118,  120,  121,  122,  123,  124,  125,  127,  128,  129,  130,  131,  132};

int exGraphData::exVertloctabRaw[245] = {   0,    2,    0,    3,    0,    3,    0,    3,    0,    3,    0,    3,    6,    0,    3,    0,    3,    0,    3,    0,    3,    0,    2,    0,    3,    7,    0,    4,    0,    4,    0,    4,    0,    4,    0,    4,    8,    0,    4,    0,    4,    0,    4,    0,    3,    0,    3,    0,    4,    8,    0,    4,    0,    4,    0,    4,    0,    4,    0,    4,    0,    4,    8,    0,    4,    0,    3,    0,    3,    0,    4,    0,    4,    8,    0,    4,    0,    4,    0,    4,    0,    4,    0,    4,    0,    4,    8,    0,    3,    0,    3,    0,    4,    0,    4,    0,    4,    0,    4,    8,    0,    4,    0,    4,    0,    4,    0,    4,    0,    4,    7,    0,    3,    0,    4,    0,    4,    0,    4,    0,    4,    0,    4,    8,    0,    4,    0,    4,    0,    4,    0,    4,    0,    3,    0,    3,    7,    0,    4,    0,    4,    0,    4,    0,    4,    0,    4,    8,    0,    4,    0,    4,    0,    4,    0,    3,    0,    3,    0,    4,    8,    0,    4,    0,    4,    0,    4,    0,    4,    0,    4,    0,    4,    8,    0,    4,    0,    3,    0,    3,    0,    4,    0,    4,    8,    0,    4,    0,    4,    0,    4,    0,    4,    0,    4,    0,    4,    8,    0,    3,    0,    3,    0,    4,    0,    4,    0,    4,    0,    4,    8,    0,    4,    0,    4,    0,    4,    0,    4,    0,    4,    7,    0,    2,    0,    3,    0,    3,    0,    3,    0,    3,    0,    3,    6,    0,    3,    0,    3,    0,    3,    0,    3,    0,    2};
int exGraphData::exVertwgtcmpRaw[245] = {1100, 1101, 1110, 1111, 1120, 1121, 1130, 1131, 1140, 1141, 1150, 1151, 1152, 1160, 1161, 1170, 1171, 1180, 1181, 1190, 1191, 1200, 1201, 1210, 1211, 1212, 1220, 1221, 1230, 1231, 1240, 1241, 1250, 1251, 1260, 1261, 1262, 1270, 1271, 1280, 1281, 1290, 1291, 1300, 1301, 1310, 1311, 1320, 1321, 1322, 1330, 1331, 1340, 1341, 1350, 1351, 1360, 1361, 1370, 1371, 1380, 1381, 1382, 1390, 1391, 1400, 1401, 1410, 1411, 1420, 1421, 1430, 1431, 1432, 1440, 1441, 1450, 1451, 1460, 1461, 1470, 1471, 1480, 1481, 1490, 1491, 1492, 1500, 1501, 1510, 1511, 1520, 1521, 1530, 1531, 1540, 1541, 1550, 1551, 1552, 1560, 1561, 1570, 1571, 1580, 1581, 1590, 1591, 1600, 1601, 1602, 1610, 1611, 1620, 1621, 1630, 1631, 1640, 1641, 1650, 1651, 1660, 1661, 1662, 1670, 1671, 1680, 1681, 1690, 1691, 1700, 1701, 1710, 1711, 1720, 1721, 1722, 1730, 1731, 1740, 1741, 1750, 1751, 1760, 1761, 1770, 1771, 1772, 1780, 1781, 1790, 1791, 1800, 1801, 1810, 1811, 1820, 1821, 1830, 1831, 1832, 1840, 1841, 1850, 1851, 1860, 1861, 1870, 1871, 1880, 1881, 1890, 1891, 1892, 1900, 1901, 1910, 1911, 1920, 1921, 1930, 1931, 1940, 1941, 1942, 1950, 1951, 1960, 1961, 1970, 1971, 1980, 1981, 1990, 1991, 2000, 2001, 2002, 2010, 2011, 2020, 2021, 2030, 2031, 2040, 2041, 2050, 2051, 2060, 2061, 2062, 2070, 2071, 2080, 2081, 2090, 2091, 2100, 2101, 2110, 2111, 2112, 2120, 2121, 2130, 2131, 2140, 2141, 2150, 2151, 2160, 2161, 2170, 2171, 2172, 2180, 2181, 2190, 2191, 2200, 2201, 2210, 2211, 2220, 2221};
int exGraphData::exVertwgtvisRaw[245] = {1200, 1201, 1210, 1211, 1220, 1221, 1230, 1231, 1240, 1241, 1250, 1251, 1252, 1260, 1261, 1270, 1271, 1280, 1281, 1290, 1291, 1300, 1301, 1310, 1311, 1312, 1320, 1321, 1330, 1331, 1340, 1341, 1350, 1351, 1360, 1361, 1362, 1370, 1371, 1380, 1381, 1390, 1391, 1400, 1401, 1410, 1411, 1420, 1421, 1422, 1430, 1431, 1440, 1441, 1450, 1451, 1460, 1461, 1470, 1471, 1480, 1481, 1482, 1490, 1491, 1500, 1501, 1510, 1511, 1520, 1521, 1530, 1531, 1532, 1540, 1541, 1550, 1551, 1560, 1561, 1570, 1571, 1580, 1581, 1590, 1591, 1592, 1600, 1601, 1610, 1611, 1620, 1621, 1630, 1631, 1640, 1641, 1650, 1651, 1652, 1660, 1661, 1670, 1671, 1680, 1681, 1690, 1691, 1700, 1701, 1702, 1710, 1711, 1720, 1721, 1730, 1731, 1740, 1741, 1750, 1751, 1760, 1761, 1762, 1770, 1771, 1780, 1781, 1790, 1791, 1800, 1801, 1810, 1811, 1820, 1821, 1822, 1830, 1831, 1840, 1841, 1850, 1851, 1860, 1861, 1870, 1871, 1872, 1880, 1881, 1890, 1891, 1900, 1901, 1910, 1911, 1920, 1921, 1930, 1931, 1932, 1940, 1941, 1950, 1951, 1960, 1961, 1970, 1971, 1980, 1981, 1990, 1991, 1992, 2000, 2001, 2010, 2011, 2020, 2021, 2030, 2031, 2040, 2041, 2042, 2050, 2051, 2060, 2061, 2070, 2071, 2080, 2081, 2090, 2091, 2100, 2101, 2102, 2110, 2111, 2120, 2121, 2130, 2131, 2140, 2141, 2150, 2151, 2160, 2161, 2162, 2170, 2171, 2180, 2181, 2190, 2191, 2200, 2201, 2210, 2211, 2212, 2220, 2221, 2230, 2231, 2240, 2241, 2250, 2251, 2260, 2261, 2270, 2271, 2272, 2280, 2281, 2290, 2291, 2300, 2301, 2310, 2311, 2320, 2321};

int exGraphData::exEdgeloctabRaw[482] = {   1,   12,    0,    2,   13,    1,    3,   14,    2,    4,   15,    3,    5,   16,    4,    6,   17,    5,    7,   18,    6,    8,   19,    7,    9,   20,    8,   10,   21,    9,   11,   22,   10,   23,    0,   13,   24,    1,   12,   14,   25,    2,   13,   15,   26,    3,   14,   16,   27,    4,   15,   17,   28,    5,   16,   18,   29,    6,   17,   19,   30,    7,   18,   20,   31,    8,   19,   21,   32,    9,   20,   22,   33,   10,   21,   23,   34,   11,   22,   35,   12,   25,   36,   13,   24,   26,   37,   14,   25,   27,   38,   15,   26,   28,   39,   16,   27,   29,   40,   17,   28,   30,   41,   18,   29,   31,   42,   19,   30,   32,   43,   20,   31,   33,   44,   21,   32,   34,   45,   22,   33,   35,   46,   23,   34,   47,   24,   37,   48,   25,   36,   38,   49,   26,   37,   39,   50,   27,   38,   40,   51,   28,   39,   41,   52,   29,   40,   42,   53,   30,   41,   43,   54,   31,   42,   44,   55,   32,   43,   45,   56,   33,   44,   46,   57,   34,   45,   47,   58,   35,   46,   59,   36,   49,   60,   37,   48,   50,   61,   38,   49,   51,   62,   39,   50,   52,   63,   40,   51,   53,   64,   41,   52,   54,   65,   42,   53,   55,   66,   43,   54,   56,   67,   44,   55,   57,   68,   45,   56,   58,   69,   46,   57,   59,   70,   47,   58,   71,   48,   61,   72,   49,   60,   62,   73,   50,   61,   63,   74,   51,   62,   64,   75,   52,   63,   65,   76,   53,   64,   66,   77,   54,   65,   67,   78,   55,   66,   68,   79,   56,   67,   69,   80,   57,   68,   70,   81,   58,   69,   71,   82,   59,   70,   83,   60,   73,   84,   61,   72,   74,   85,   62,   73,   75,   86,   63,   74,   76,   87,   64,   75,   77,   88,   65,   76,   78,   89,   66,   77,   79,   90,   67,   78,   80,   91,   68,   79,   81,   92,   69,   80,   82,   93,   70,   81,   83,   94,   71,   82,   95,   72,   85,   96,   73,   84,   86,   97,   74,   85,   87,   98,   75,   86,   88,   99,   76,   87,   89,  100,   77,   88,   90,  101,   78,   89,   91,  102,   79,   90,   92,  103,   80,   91,   93,  104,   81,   92,   94,  105,   82,   93,   95,  106,   83,   94,  107,   84,   97,  108,   85,   96,   98,  109,   86,   97,   99,  110,   87,   98,  100,  111,   88,   99,  101,  112,   89,  100,  102,  113,   90,  101,  103,  114,   91,  102,  104,  115,   92,  103,  105,  116,   93,  104,  106,  117,   94,  105,  107,  118,   95,  106,  119,   96,  109,  120,   97,  108,  110,  121,   98,  109,  111,  122,   99,  110,  112,  123,  100,  111,  113,  124,  101,  112,  114,  125,  102,  113,  115,  126,  103,  114,  116,  127,  104,  115,  117,  128,  105,  116,  118,  129,  106,  117,  119,  130,  107,  118,  131,  108,  121,  109,  120,  122,  110,  121,  123,  111,  122,  124,  112,  123,  125,  113,  124,  126,  114,  125,  127,  115,  126,  128,  116,  127,  129,  117,  128,  130,  118,  129,  131,  119,  130};
int exGraphData::exEdgewgtcmpRaw[482] = {1100, 1101, 1100, 1102, 1103, 1102, 1104, 1105, 1104, 1106, 1107, 1106, 1108, 1109, 1108, 1110, 1111, 1110, 1112, 1113, 1112, 1114, 1115, 1114, 1116, 1117, 1116, 1118, 1119, 1118, 1120, 1121, 1120, 1122, 1101, 1123, 1124, 1103, 1123, 1125, 1126, 1105, 1125, 1127, 1128, 1107, 1127, 1129, 1130, 1109, 1129, 1131, 1132, 1111, 1131, 1133, 1134, 1113, 1133, 1135, 1136, 1115, 1135, 1137, 1138, 1117, 1137, 1139, 1140, 1119, 1139, 1141, 1142, 1121, 1141, 1143, 1144, 1122, 1143, 1145, 1124, 1146, 1147, 1126, 1146, 1148, 1149, 1128, 1148, 1150, 1151, 1130, 1150, 1152, 1153, 1132, 1152, 1154, 1155, 1134, 1154, 1156, 1157, 1136, 1156, 1158, 1159, 1138, 1158, 1160, 1161, 1140, 1160, 1162, 1163, 1142, 1162, 1164, 1165, 1144, 1164, 1166, 1167, 1145, 1166, 1168, 1147, 1169, 1170, 1149, 1169, 1171, 1172, 1151, 1171, 1173, 1174, 1153, 1173, 1175, 1176, 1155, 1175, 1177, 1178, 1157, 1177, 1179, 1180, 1159, 1179, 1181, 1182, 1161, 1181, 1183, 1184, 1163, 1183, 1185, 1186, 1165, 1185, 1187, 1188, 1167, 1187, 1189, 1190, 1168, 1189, 1191, 1170, 1192, 1193, 1172, 1192, 1194, 1195, 1174, 1194, 1196, 1197, 1176, 1196, 1198, 1199, 1178, 1198, 1200, 1201, 1180, 1200, 1202, 1203, 1182, 1202, 1204, 1205, 1184, 1204, 1206, 1207, 1186, 1206, 1208, 1209, 1188, 1208, 1210, 1211, 1190, 1210, 1212, 1213, 1191, 1212, 1214, 1193, 1215, 1216, 1195, 1215, 1217, 1218, 1197, 1217, 1219, 1220, 1199, 1219, 1221, 1222, 1201, 1221, 1223, 1224, 1203, 1223, 1225, 1226, 1205, 1225, 1227, 1228, 1207, 1227, 1229, 1230, 1209, 1229, 1231, 1232, 1211, 1231, 1233, 1234, 1213, 1233, 1235, 1236, 1214, 1235, 1237, 1216, 1238, 1239, 1218, 1238, 1240, 1241, 1220, 1240, 1242, 1243, 1222, 1242, 1244, 1245, 1224, 1244, 1246, 1247, 1226, 1246, 1248, 1249, 1228, 1248, 1250, 1251, 1230, 1250, 1252, 1253, 1232, 1252, 1254, 1255, 1234, 1254, 1256, 1257, 1236, 1256, 1258, 1259, 1237, 1258, 1260, 1239, 1261, 1262, 1241, 1261, 1263, 1264, 1243, 1263, 1265, 1266, 1245, 1265, 1267, 1268, 1247, 1267, 1269, 1270, 1249, 1269, 1271, 1272, 1251, 1271, 1273, 1274, 1253, 1273, 1275, 1276, 1255, 1275, 1277, 1278, 1257, 1277, 1279, 1280, 1259, 1279, 1281, 1282, 1260, 1281, 1283, 1262, 1284, 1285, 1264, 1284, 1286, 1287, 1266, 1286, 1288, 1289, 1268, 1288, 1290, 1291, 1270, 1290, 1292, 1293, 1272, 1292, 1294, 1295, 1274, 1294, 1296, 1297, 1276, 1296, 1298, 1299, 1278, 1298, 1300, 1301, 1280, 1300, 1302, 1303, 1282, 1302, 1304, 1305, 1283, 1304, 1306, 1285, 1307, 1308, 1287, 1307, 1309, 1310, 1289, 1309, 1311, 1312, 1291, 1311, 1313, 1314, 1293, 1313, 1315, 1316, 1295, 1315, 1317, 1318, 1297, 1317, 1319, 1320, 1299, 1319, 1321, 1322, 1301, 1321, 1323, 1324, 1303, 1323, 1325, 1326, 1305, 1325, 1327, 1328, 1306, 1327, 1329, 1308, 1330, 1310, 1330, 1331, 1312, 1331, 1332, 1314, 1332, 1333, 1316, 1333, 1334, 1318, 1334, 1335, 1320, 1335, 1336, 1322, 1336, 1337, 1324, 1337, 1338, 1326, 1338, 1339, 1328, 1339, 1340, 1329, 1340};
int exGraphData::exEdgewgtvisRaw[482] = {1200, 1201, 1200, 1202, 1203, 1202, 1204, 1205, 1204, 1206, 1207, 1206, 1208, 1209, 1208, 1210, 1211, 1210, 1212, 1213, 1212, 1214, 1215, 1214, 1216, 1217, 1216, 1218, 1219, 1218, 1220, 1221, 1220, 1222, 1201, 1223, 1224, 1203, 1223, 1225, 1226, 1205, 1225, 1227, 1228, 1207, 1227, 1229, 1230, 1209, 1229, 1231, 1232, 1211, 1231, 1233, 1234, 1213, 1233, 1235, 1236, 1215, 1235, 1237, 1238, 1217, 1237, 1239, 1240, 1219, 1239, 1241, 1242, 1221, 1241, 1243, 1244, 1222, 1243, 1245, 1224, 1246, 1247, 1226, 1246, 1248, 1249, 1228, 1248, 1250, 1251, 1230, 1250, 1252, 1253, 1232, 1252, 1254, 1255, 1234, 1254, 1256, 1257, 1236, 1256, 1258, 1259, 1238, 1258, 1260, 1261, 1240, 1260, 1262, 1263, 1242, 1262, 1264, 1265, 1244, 1264, 1266, 1267, 1245, 1266, 1268, 1247, 1269, 1270, 1249, 1269, 1271, 1272, 1251, 1271, 1273, 1274, 1253, 1273, 1275, 1276, 1255, 1275, 1277, 1278, 1257, 1277, 1279, 1280, 1259, 1279, 1281, 1282, 1261, 1281, 1283, 1284, 1263, 1283, 1285, 1286, 1265, 1285, 1287, 1288, 1267, 1287, 1289, 1290, 1268, 1289, 1291, 1270, 1292, 1293, 1272, 1292, 1294, 1295, 1274, 1294, 1296, 1297, 1276, 1296, 1298, 1299, 1278, 1298, 1300, 1301, 1280, 1300, 1302, 1303, 1282, 1302, 1304, 1305, 1284, 1304, 1306, 1307, 1286, 1306, 1308, 1309, 1288, 1308, 1310, 1311, 1290, 1310, 1312, 1313, 1291, 1312, 1314, 1293, 1315, 1316, 1295, 1315, 1317, 1318, 1297, 1317, 1319, 1320, 1299, 1319, 1321, 1322, 1301, 1321, 1323, 1324, 1303, 1323, 1325, 1326, 1305, 1325, 1327, 1328, 1307, 1327, 1329, 1330, 1309, 1329, 1331, 1332, 1311, 1331, 1333, 1334, 1313, 1333, 1335, 1336, 1314, 1335, 1337, 1316, 1338, 1339, 1318, 1338, 1340, 1341, 1320, 1340, 1342, 1343, 1322, 1342, 1344, 1345, 1324, 1344, 1346, 1347, 1326, 1346, 1348, 1349, 1328, 1348, 1350, 1351, 1330, 1350, 1352, 1353, 1332, 1352, 1354, 1355, 1334, 1354, 1356, 1357, 1336, 1356, 1358, 1359, 1337, 1358, 1360, 1339, 1361, 1362, 1341, 1361, 1363, 1364, 1343, 1363, 1365, 1366, 1345, 1365, 1367, 1368, 1347, 1367, 1369, 1370, 1349, 1369, 1371, 1372, 1351, 1371, 1373, 1374, 1353, 1373, 1375, 1376, 1355, 1375, 1377, 1378, 1357, 1377, 1379, 1380, 1359, 1379, 1381, 1382, 1360, 1381, 1383, 1362, 1384, 1385, 1364, 1384, 1386, 1387, 1366, 1386, 1388, 1389, 1368, 1388, 1390, 1391, 1370, 1390, 1392, 1393, 1372, 1392, 1394, 1395, 1374, 1394, 1396, 1397, 1376, 1396, 1398, 1399, 1378, 1398, 1400, 1401, 1380, 1400, 1402, 1403, 1382, 1402, 1404, 1405, 1383, 1404, 1406, 1385, 1407, 1408, 1387, 1407, 1409, 1410, 1389, 1409, 1411, 1412, 1391, 1411, 1413, 1414, 1393, 1413, 1415, 1416, 1395, 1415, 1417, 1418, 1397, 1417, 1419, 1420, 1399, 1419, 1421, 1422, 1401, 1421, 1423, 1424, 1403, 1423, 1425, 1426, 1405, 1425, 1427, 1428, 1406, 1427, 1429, 1408, 1430, 1410, 1430, 1431, 1412, 1431, 1432, 1414, 1432, 1433, 1416, 1433, 1434, 1418, 1434, 1435, 1420, 1435, 1436, 1422, 1436, 1437, 1424, 1437, 1438, 1426, 1438, 1439, 1428, 1439, 1440, 1429, 1440};
