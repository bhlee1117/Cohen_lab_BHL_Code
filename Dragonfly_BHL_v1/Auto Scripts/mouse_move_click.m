function mouse_move_click(xy)
import java.awt.Robot
import java.awt.event.InputEvent
input = Robot;
leftbutton = java.awt.event.InputEvent.BUTTON1_MASK;

input.mouseMove(0, 0);
input.mouseMove(xy(1), xy(2));
pause(.1)
input.mousePress(leftbutton);
input.mouseRelease(leftbutton);