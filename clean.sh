find ./build ! -name .gitkeep ! -wholename ./build -delete
rm -rf Debug
rm -rf Release
find . -wholename ./spirit* -delete
find ./core/python/Spirit -mindepth 1 -name *Spirit* -delete
find ./core/julia/Spirit  -mindepth 1 -name *Spirit* -delete
find ./ui-web  -mindepth 1 -name spirit.js* -delete
