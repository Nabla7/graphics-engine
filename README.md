## Gequoteerde functionaliteit

V: Werkend  
-: Deels werkend met gekende problemen (onderaan beschreven)  
X: Niet werkend of niet geïmplementeerd  
(blanco): TODO  


|   | Functionaliteit      | Status |
|---|---------------------------|---|
| 1 | 2D L-systemen             | V |
|   | Met haakjes               | V |
|   | Stochastisch              | X |
| 2 | Transformaties            | V |
|   | Eye-point                 | V |
|   | Projectie                 | V |
| 3 | Platonische Lichamen      | V |
|   | Kegel en cylinder         | V |
|   | Bol                       | V |
|   | Torus                     | V |
|   | 3D L-systemen             | X |
| 4 | Z-buffering (lijnen)      | - |
| 5 | Triangulatie              | V |
|   | Z-buffering (driehoeken)  | - |
| 6 | 3D fractalen              | V |
|   | BuckyBall                 | X |
|   | Mengerspons               | X |
|   | View Frustum              | X |
| 7 | Ambient licht             | X |
|   | Diffuus licht (oneindig)  | X |
|   | Diffuus licht (puntbron)  | X |
|   | Speculair licht           | X |
| 8 | Schaduw                   | X |
|   | Texture mapping           | X |
| 9 | Bollen en cylinders       | X |
|   | UV-coordinaten            | X |
|   | Cube mapping              | X |

Geïmplementeerde vorm van texture mapping: ...

## Gekende problemen 

Voor z-buffering (lijnen) is de output duidelijk beter dan de output zonder z-buffering maar voor objecten met veel lijnen is de performance niet even goed als de output van de reference implementatie.

Voor Z-buffering (driehoeken) werkt alles even goed als de reference implementatie behalve voor de z_buffering079.ini test case


## Niet-gequoteerde functionaliteit
...

## Extra functionaliteit, niet in de opgaves beschreven
...

