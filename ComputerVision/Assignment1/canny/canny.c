#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAX_VALUE 255
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

struct matrix{
    int** data;
    int width;
    int height;
};

struct matrix zeroedMatrix(int width, int height){
    struct matrix mat;
    
    mat.width = width;
    mat.height = height;
    mat.data = calloc(height, sizeof(int*));
    for(int i = 0; i < mat.height; i++)
        mat.data[i] = calloc(width, sizeof(int));
    
    return mat;
}

void freeMatrix(struct matrix mat){
    for(int i = 0; i < mat.height; i++)
        free(mat.data[i]);
    free(mat.data);
}

int floatToInt(double d){
    return (int)floor(d + 0.5 - 1e-9);
}

struct matrix normalize(struct matrix mat){
    struct matrix normal = zeroedMatrix(mat.width, mat.height);
    
    int min = 1e9;
    int max = 0;
    for(int i = 0; i < mat.height; i++){
        for(int j = 0; j < mat.width; j++){
            max = MAX(max, mat.data[i][j]);
            min = MIN(min, mat.data[i][j]);
        }
    }
    
    int range = max - min;
    double factor = MAX_VALUE / (double)range;
    for(int i = 0; i < mat.height; i++)
        for(int j = 0; j < mat.width; j++)
            normal.data[i][j] = MIN(floatToInt(factor * (mat.data[i][j] - min)), MAX_VALUE);
        
    return normal;
}

void printMatrix(struct matrix mat){
    for(int i = 0; i < mat.height; i++){
        for(int j = 0; j < mat.width; j++)
            printf("%d ", mat.data[i][j]);
        printf("\n");
    }
}

struct matrix readPgm(char* pgmFileName){
    
    struct matrix image;
    
    char line[100];
    FILE* fp = fopen(pgmFileName, "rb");
    
    fgets(line, 100, fp);
    if(strlen(line) == 3)
        fgets(line, 100, fp);
    
    int width, height;
    sscanf(line, "%d %d", &width, &height);
    image = zeroedMatrix(width, height);
    
    int throwaway;
    fgets(line, 100, fp);
    sscanf(line, "%d", &throwaway);
    
    for(int i = 0; i < image.height; i++)
        for(int j = 0; j < image.width; j++)
            image.data[i][j] = fgetc(fp);
    
    fclose(fp);
    return image;
}

struct matrix getGaussianFilter(int horizontal, double sigma){
    int center = 3 * floatToInt(sigma);
    int dim = 2 * center + 1;
    double factor = 1 / (center * exp((-center * center) / (sigma * sigma)));
    
    struct matrix filter = zeroedMatrix(dim, dim);
        
    for(int i = 0; i < filter.height; i++){
        int y = i - center;
        for(int j = 0; j < filter.width; j++){
            int x = j - center;
            double gauss = exp(-(x * x + y * y)/(2 * sigma * sigma));
            if(horizontal)
                filter.data[i][j] = floatToInt(x * gauss * factor);
            else
                filter.data[i][j] = floatToInt(y * gauss * factor);
        }
    }
    
    return filter;
}

int weightedSum(struct matrix image, struct matrix filter, int row, int col){
    int result = 0;
    for(int i = 0; i < filter.height; i++)
        for(int j = 0; j < filter.width; j++)
            result += image.data[row + i][col + j] * filter.data[i][j];
    return result;
}

struct matrix applyConvolution(struct matrix image, struct matrix filter){
    struct matrix convolution = zeroedMatrix(image.width, image.height);
    struct matrix paddedImage = zeroedMatrix(image.width + filter.width - 1, image.height + filter.height - 1);
    
    for(int i = 0; i < image.height; i++)
        for(int j = 0; j < image.width; j++)
            paddedImage.data[i + (filter.height - 1) / 2][j + (filter.width - 1) / 2] = image.data[i][j];
    for(int i = 0; i < convolution.height; i++)
        for(int j = 0; j < convolution.width; j++)
            convolution.data[i][j] = weightedSum(paddedImage, filter, i, j);
    
    freeMatrix(paddedImage);
    
    return convolution;
}

int magnitude(int x, int y){
    return floatToInt(sqrt(x * x + y * y));
}

struct matrix buildMagnitudeMap(struct matrix h_conv, struct matrix v_conv){
    struct matrix magMap = zeroedMatrix(h_conv.width, h_conv.height);
    
    for(int i = 0; i < magMap.height; i++)
        for(int j = 0; j < magMap.width; j++)
            magMap.data[i][j] = magnitude(h_conv.data[i][j], v_conv.data[i][j]);
    return magMap;
}

struct matrix extractPeakMap(struct matrix magMap){
    struct matrix peakMap = zeroedMatrix(magMap.width, magMap, height);
    int offset[4][2] = {
        {0, 1},
        {0, -1},
        {1, 0},
        {-1, 0}
    };
    
    for(int i = 1; i < peakMap.height - 1; i++){
        for(int j = 1; j < peakMap.width - 1; j++){
            int isPeak = 1;
            for(int k = 0; k < 4; k++){
                if(magMap.data[i][j] < magMap.data[i + offset[k][0]][j + offset[k][1]]){
                    isPeak = 0;
                    break;
                }
            }
            peakMap.data[i][j] = isPeak;
        }
    }
    
    return peakMap;
}

void savePgm(char* imageFile, struct matrix image){
    image = normalize(image);
    
    FILE* fp = fopen(imageFile, "wb");
    fprintf(fp, "P5\n");
    fprintf(fp, "%d %d\n", image.width, image.height);
    fprintf(fp, "%d\n", MAX_VALUE);
    for(int i = 0; i < image.height; i++)
        for(int j = 0; j < image.width; j++)
            fprintf(fp, "%c", image.data[i][j]);
    
    fclose(fp);
}

int main(int argc, char** argv){
    char* imageName = argv[1];
    double sigma = atof(argv[2]);
    
    printf("Reading Image Data...\n");
    struct matrix image = readPgm(imageName);
    
    printf("Applying Convolutions...\n");
    struct matrix h_conv = applyConvolution(image, getGaussianFilter(1, sigma));
    struct matrix v_conv = applyConvolution(image, getGaussianFilter(0, sigma));
    
    printf("Building Magnitude Map...\n");
    struct matrix magMap = buildMagnitudeMap(h_conv, v_conv);
    
    printf("Saving Images...\n");
    
    savePgm("horizontalGrad.pgm", h_conv);
    savePgm("verticalGrad.pgm", v_conv);
    savePgm("magnitudeMap.pgm", magMap);
    
    freeMatrix(h_conv);
    freeMatrix(v_conv);
    freeMatrix(magMap);
    
    return 0;
}