## Ускорение видео

```bash
ffmpeg -i input.mp4 -filter:v "setpts=0.5*PTS" output.mp4
```