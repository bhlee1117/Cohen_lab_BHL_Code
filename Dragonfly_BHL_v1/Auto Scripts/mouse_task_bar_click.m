function mouse_task_bar_click(xy)
import java.awt.Robot
import java.awt.event.InputEvent
input = Robot;
leftbutton = java.awt.event.InputEvent.BUTTON1_MASK;

xy_tb = xy;
xy_window = xy(1)+[70 -70];

input.mouseMove(0, 0);
input.mouseMove(xy_tb(1), xy_tb(2));
pause(.2)
input.mousePress(leftbutton);
input.mouseRelease(leftbutton);

pause(.1)
input.mouseMove(0, 0);
input.mouseMove(xy_window(1), xy_window(2));
input.mousePress(leftbutton);
input.mouseRelease(leftbutton);