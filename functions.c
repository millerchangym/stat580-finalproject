int stringNumericCheck (char *str) {
        int i = 0;
        int notNumericCount = 0;

        for (i = 0; i < strlen(str); i++) {
                notNumericCount += !(isdigit(str[i]));
        }

        return !(notNumericCount > 0);
}
