//
//  AppDelegate.h
//  StatCalc
//
//  Created by Koichi Oshio on 12/25/13.
//
//

#import <Cocoa/Cocoa.h>
#import <NumKit/NumKit.h>

@interface AppDelegate : NSObject <NSApplicationDelegate>
{
// t-table
    IBOutlet NSTextField    *t_dfField;
    IBOutlet NSTextField    *t_Field;
    IBOutlet NSTextField    *t_pField;
    IBOutlet NSTextField    *t_p05Field;
    IBOutlet NSTextField    *t_p01Field;
// F-table
    IBOutlet NSTextField    *F_df1Field;
    IBOutlet NSTextField    *F_df2Field;
    IBOutlet NSTextField    *F_Field;
    IBOutlet NSTextField    *F_pField;
    IBOutlet NSTextField    *F_p05Field;
    IBOutlet NSTextField    *F_p01Field;
    IBOutlet NSTextField    *F_p05LField;
    IBOutlet NSTextField    *F_p01LField;
}

@property (assign) IBOutlet NSWindow *window;

// t-table
- (IBAction)updateT:(id)sender;

// F-table
- (IBAction)updateF:(id)sender;

@end
