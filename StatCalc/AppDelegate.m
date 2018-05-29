//
//  AppDelegate.m
//  StatCalc
//
//  Created by Koichi Oshio on 12/25/13.
//
//

#import "AppDelegate.h"

@implementation AppDelegate

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification
{
    NSNumberFormatter   *format1 = [[NSNumberFormatter alloc] init];
    NSNumberFormatter   *format2 = [[NSNumberFormatter alloc] init];

    [format1    setFormat:@"#0.####"];
    [t_p05Field setFormatter:format1];
    [t_p01Field setFormatter:format1];
    [F_p05Field setFormatter:format1];
    [F_p01Field setFormatter:format1];
    [F_p05LField setFormatter:format1];
    [F_p01LField setFormatter:format1];

    [format2    setFormat:@"#0.####"];
    [t_pField   setFormatter:format2];
    [F_pField   setFormatter:format2];

    [self updateT:self];
    [self updateF:self];
}

// t-table
- (IBAction)updateT:(id)sender
{
    float   t = [t_Field floatValue];
    int     df = [t_dfField intValue];
    float   p, p01, p05;

    if (df < 1) {
        df = 1;
        [t_dfField setIntValue:df];
    }
    p = Num_t2p(t, df);
    [t_pField setFloatValue:p]; // one-sided

    p05 = Num_tlimit(0.05, df);    // Num_tlimit returns two-sided value
    [t_p05Field setFloatValue:p05];

    p01 = Num_tlimit(0.01, df);    //  Num_tlimit returns two-sided value
    [t_p01Field setFloatValue:p01];
}

// F-table
- (IBAction)updateF:(id)sender
{
    float   F = [F_Field floatValue];
    int     df1 = [F_df1Field intValue];
    int     df2 = [F_df2Field intValue];
    float   p;

    if (F <= 0) {
        F = 1.0;
        [F_Field setFloatValue:F];
    }
    if (df1 < 2) {
        df1 = 2;
        [F_df1Field setIntValue:df1];
    }
    if (df2 < 2) {
        df2 = 2;
        [F_df2Field setIntValue:df2];
    }
    p = Num_f2p(F, df1, df2);
    [F_pField setFloatValue:p];

    p = Num_flimit(0.05, df1, df2, NO);
    [F_p05Field setFloatValue:p];
    p = Num_flimit(0.05, df1, df2, YES);
    [F_p05LField setFloatValue:p];

    p = Num_flimit(0.01, df1, df2, NO);
    [F_p01Field setFloatValue:p];
    p = Num_flimit(0.01, df1, df2, YES);
    [F_p01LField setFloatValue:p];
}

@end
