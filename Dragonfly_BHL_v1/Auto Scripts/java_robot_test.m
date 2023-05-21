import java.awt.Robot
import java.awt.event.InputEvent
import java.awt.event.KeyEvent
input = Robot;
leftbutton = java.awt.event.InputEvent.BUTTON1_MASK;
kctrl = java.awt.event.KeyEvent.VK_CONTROL;
kv = java.awt.event.KeyEvent.VK_V;
kshift = java.awt.event.KeyEvent.VK_SHIFT;
kbs = java.awt.event.KeyEvent.VK_BACK_SPACE;
%% spinview taskbar location
input.mouseMove(0, 0);
input.mouseMove(1000, 1170);
pause(.2)
input.mousePress(leftbutton);
input.mouseRelease(leftbutton);
%%
pause(.1)
input.mouseMove(0, 0);
input.mouseMove(1020, 1100);
input.mousePress(leftbutton);
input.mouseRelease(leftbutton);

%% start rec location
input.mouseMove(0, 0);
input.mouseMove(1850, 730);
%% stop rec location
input.mouseMove(0, 0);
input.mouseMove(1750, 730);
% input.mousePress(InputEvent.BUTTON1_MASK);
% input.mouseRelease(InputEvent.BUTTON1_MASK);
%% text area location
input.mouseMove(0, 0);
input.mouseMove(1700, 700);
input.mouseMove(1548, 80);

input.mousePress(leftbutton);
input.keyPress(kshift);
input.mouseMove(1850, 150);
input.mouseRelease(leftbutton);
input.keyRelease(kshift);
%%
input.keyPress(kbs)
input.keyRelease(kbs)

clipboard('copy','E:\M')
input.keyPress(kctrl)
input.keyPress(kv)
input.keyRelease(kv)
input.keyRelease(kctrl)
%%
input.keyRelease(kbs)
input.keyRelease(kbs)
%%
input.mousePress(leftbutton);
input.mouseRelease(leftbutton);
%%
input.mousePress(java.awt.event.InputEvent.BUTTON1_MASK);