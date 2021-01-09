//
//  NumMetal.m
//  NumKit
//
//  Created by Koichi Oshio on 2021/01/02.
//

#import "NumKit.h"

@implementation NumMetalDevice

- (id) initWithLibrary:(NSString *)path
{
    NSError     *err;
    NSBundle    *bundle;

    self = [super init];
    if (self) {
        device = MTLCreateSystemDefaultDevice();
        if (device == nil) return nil;
 
        if (path) {
            library = [device newLibraryWithFile:path error:&err];
        } else {
            // bundle or URL works, path doesn't work
            bundle = [NSBundle bundleForClass:[self class]];
            library = [device newDefaultLibraryWithBundle:bundle error:&err];
        }
        if (library == nil) {
            printf("library == nil\n");
            exit(0);
        }
    }
    queue = [device newCommandQueue];   // command buffer needs to be created every time

    return self;
}    

- (void)setFunction:(NSString *)funcName;
{
    NSError *err;
    function = [library newFunctionWithName:funcName];
    pipeline = [device newComputePipelineStateWithFunction:function error:&err];

}

- (id)encoder
{
    return encoder;
}

- (id<MTLBuffer>) bufferWithLength:(NSUInteger)len options:(MTLResourceOptions)opt
{
    return [device newBufferWithLength:len options:opt];
}

- (void)setBuffer:(id<MTLBuffer>)buffer offset:(NSInteger)ofs atIndex:(NSInteger)ix
{
    [encoder setBuffer:buffer offset:ofs atIndex:ix];
}

- (void)commit
{
    [commandBuffer commit];
}

- (void)waitUntilCompleted
{
    [commandBuffer waitUntilCompleted];
}

- (NSUInteger)maxTotalThreadsPerThreadgroup
{
    return [pipeline maxTotalThreadsPerThreadgroup];
}

- (void)dispatchThreads:(MTLSize)gridSize threadsPerThreadgroup:(MTLSize)threadgroupSize
{
    [encoder dispatchThreads:gridSize threadsPerThreadgroup:threadgroupSize];
}

- (void)startEncoding
{
    commandBuffer = [queue commandBuffer];
    encoder = [commandBuffer computeCommandEncoder];
    [encoder setComputePipelineState:pipeline];
}

- (void)endEncoding
{
    [encoder endEncoding];
}

+ (void)query
{
    id <MTLDevice>      device;
    MTLSize             size;

    device = MTLCreateSystemDefaultDevice();
    if (device == nil) return;

    printf("maxBufferLength = %ld\n", [device maxBufferLength]);
    printf("recommendedMaxWorkingSetSize = %lld\n", [device recommendedMaxWorkingSetSize]);
    size = [device maxThreadsPerThreadgroup];
    printf("maxThreadsPerThreadGroup (%ld, %ld, %ld)\n", size.width, size.height, size.depth);
}

@end

