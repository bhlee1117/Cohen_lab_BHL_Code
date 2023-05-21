function paste_to_textarea(xy,txt)

import java.awt.Robot
import java.awt.event.InputEvent
import java.awt.event.KeyEvent
input = Robot;
leftbutton = java.awt.event.InputEvent.BUTTON1_MASK;
kshift = java.awt.event.KeyEvent.VK_SHIFT;
kctrl = java.awt.event.KeyEvent.VK_CONTROL;
kv = java.awt.event.KeyEvent.VK_V;
kbs = java.awt.event.KeyEvent.VK_BACK_SPACE;

input.mouseMove(0, 0);
input.mouseMove(1700, 700);
input.mouseMove(xy(1), xy(2));
pause(.1)
input.mousePress(leftbutton);
pause(.1)
input.keyPress(kshift)
input.mouseMove(xy(1)+350, xy(2)+80);
input.mouseRelease(leftbutton);
input.keyRelease(kshift)
pause(.1)
input.keyPress(kbs)
input.keyRelease(kbs)

clipboard('copy',txt)
pause(.1)
input.keyPress(kctrl)
input.keyPress(kv)
input.keyRelease(kv)
input.keyRelease(kctrl)